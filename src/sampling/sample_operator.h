#pragma once

#include "grid/grid_index_types.h"
#include "sampling/sample_clouds.h"
#include "grid/exchange_message.h"
#include "core/ctrs.h"

namespace spade::sampling
{
    template <typename grid_t, typename coeff_t, const grid::array_centering ctr, const int coeff_ct>
    struct interp_oper_t
    {
        using grid_type  = grid_t;
        using index_type = typename spade::grid::get_index_type<ctr>::array_type;
        using coeff_type = coeff_t;
        template <typename data_t> using arr_t = ctrs::array<data_t, coeff_ct>;
        
        static constexpr int num_coeff()                          { return coeff_ct; }
        static constexpr spade::grid::array_centering centering() { return ctr; }
        spade::device::shared_vector<arr_t<index_type>> idx;
        spade::device::shared_vector<arr_t<coeff_type>> coeff;
        
        std::size_t size() const { return idx.size(); }
        
        void transfer()
        {
            idx.transfer();
            coeff.transfer();
        }
    };
    
    template <typename grid_t, typename coeff_t, typename data_t, const grid::array_centering ctr, const int coeff_ct, typename device_t>
    struct array_interpolation_t
    {
        static constexpr int num_coeff()                          { return coeff_ct; }
        static constexpr spade::grid::array_centering centering() { return ctr; }
        
        using oper_t = interp_oper_t<grid_t, coeff_t, ctr, coeff_ct>;
        
        oper_t mine;   // Note that this will simply have zero entries for points that I don't own.
        oper_t others; // --> [0, 0, 0, 0, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5] e.g. I am rank 3 of 6
        
        device::shared_vector<std::size_t> send_sizes;
        device::shared_vector<std::size_t> scatter_offsets;
        
        // We will use aliasing for this
        grid::exchange_message_t<data_t, device_t> result_message;
        
        void transfer()
        {
            mine.transfer();
            others.transfer();
            send_sizes.transfer();
            scatter_offsets.transfer();
        }
    };
    
    template <
        typename arr_t,
        ctrs::basic_array container_t,
        typename strategy_t,
        std::invocable<typename arr_t::index_type> exclusion_crit_t>
    static auto create_interpolation(const arr_t& arr, const container_t& points, const strategy_t& strategy, const exclusion_crit_t& exclusion_crit)
    {
        static_assert(
            std::floating_point<typename arr_t::fundamental_type>,
            "cannot currently produce interpolation operator for non floating-point arrays");
        
        static_assert(
            std::same_as<typename arr_t::grid_type::coord_sys_type, coords::identity<typename arr_t::grid_type::coord_type>>,
            "interpolation not yet supported for non-cartesian coordinate systems");
        
        static_assert(
            arr_t::centering_type() == grid::cell_centered,
            "interpolation not yet supported for anything other than cell centering");
        
        using pnt_t   = typename container_t::value_type;
        using grid_t  = typename arr_t::grid_type;
        using coeff_t = typename arr_t::fundamental_type;
        
        const auto& grid  = arr.get_grid();
        const auto& group = grid.group();
        constexpr int dim = grid.dim();
        
        constexpr int num_sten_pt = strategy_t::stencil_size();
        using device_t = typename utils::remove_all<decltype(arr.device())>::type;
        using coord_float_t = typename pnt_t::value_type;
        using index_t       = typename arr_t::index_type;
        using val_t         = typename arr_t::fundamental_type;
        using output_t = array_interpolation_t<grid_t, coeff_t, val_t, arr_t::centering_type(), num_sten_pt, device_t>;
        output_t output;
        
        // We will use a BVH to check which block each point is in.
        // Might want to consider pre-computing this somewhere
        geom::bvh_t<dim, typename pnt_t::value_type> block_bvh;
        geom::bvh_params_t bvhparams{3, 1000};
        
        // Perform in computational coordinates, note that we strictly
        // need a 2D BVH if the grid is 2D
        
        const auto dx_max = grid.compute_dx_max();
        const auto dx_min = grid.compute_dx_min();
        
        const auto bbx = grid.get_bounds();
        using bnd_t = bound_box_t<typename pnt_t::value_type, dim>;
        bnd_t bnd;
        for (int d = 0; d < dim; ++d)
        {
            bnd.min(d) = bbx.min(d) - arr.get_num_exchange(d)*dx_max[d];
            bnd.max(d) = bbx.max(d) + arr.get_num_exchange(d)*dx_max[d];
        }
        
        const auto el_check = [&](const std::size_t& lb_glob, const auto& bnd_in)
        {
            const auto lb = utils::tag[partition::global](lb_glob);
            auto block_box = grid.get_bounding_box(lb);
            const auto dx = grid.get_dx(lb);
            for (int d = 0; d < dim; ++d)
            {
                block_box.min(d) = block_box.min(d) - arr.get_num_exchange(d)*dx[d];
                block_box.max(d) = block_box.max(d) + arr.get_num_exchange(d)*dx[d];
            }
            return block_box.intersects(bnd_in);
        };
        
        block_bvh.calculate(bnd, grid.get_num_global_blocks(), el_check, bvhparams);
        // We will now compute the indices of each point in the interpolation cloud using the BVH
        
        output.mine.idx.resize(points.size());
        output.mine.coeff.resize(points.size());
        
        std::vector<pnt_t> failed_points_interp;
        std::vector<pnt_t> failed_points_sampling;
        
        using cpu_t = decltype(device::cpu);
        grid::exchange_message_t<coord_float_t, cpu_t> point_msg;
        point_msg.set_size(group.size());
        
        struct lfinfo
        {
            int rank, lb;
            pnt_t pt;
            std::size_t gidx;
        };
        
        std::vector<lfinfo> off_infos;
        
        const auto& part = grid.get_partition();
        
        output.send_sizes.resize(group.size(), 0);
        output.result_message.set_size(group.size());
        
        const auto& get_glob_block = [&](const pnt_t& x_sample, bool& search_failed)
        {
            search_failed = false;
            auto lb = utils::tag[partition::global](-1);
            const auto eval = [&](const std::size_t& lb_cand)
            {
                auto lb_tagged = utils::tag[partition::global](lb_cand);
                const auto block_box = grid.get_bounding_box(lb_tagged);
                if (block_box.contains(x_sample)) lb.value = lb_cand;
            };
            
            using vec_t = typename decltype(block_bvh)::pnt_t;
            const auto x_arr = ctrs::to_array(x_sample);
            block_bvh.check_elements(eval, utils::bbox_around(x_arr, dx_min));
            if (lb.value < 0)
            {
                // If we get here, one of two things happened:
                // 1. The point we are asking for is genuinely outside of the domain
                // 2. The point we are asking for is only just outside of the domain, in the exchange cells
                // We will check for (2), and if it is not the case, assume (1) and throw an exception
                int num_checked = 0;
                const auto eval_extended = [&](const std::size_t& lb_cand)
                {
                    num_checked++;
                    auto lb_tagged = utils::tag[partition::global](lb_cand);
                    auto block_box = grid.get_bounding_box(lb_tagged);
                    const auto dx  = grid.get_dx(lb_tagged);
                    for (int d = 0; d < dx.size(); ++d)
                    {
                        block_box.min(d) -= dx[d]*arr.get_num_exchange(d);
                        block_box.max(d) += dx[d]*arr.get_num_exchange(d);
                    }
                    if (block_box.contains(x_sample)) lb.value = lb_cand;
                };
            
                using vec_t = typename decltype(block_bvh)::pnt_t;
                block_bvh.check_elements(eval_extended, ctrs::to_array(x_sample));
                if (lb.value < 0)
                {
                    search_failed = true;
                }
            }
            
            return lb;
        };
        
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            bool failed = false;
            const auto& x_sample = points[i];
            const auto lb = get_glob_block(x_sample, failed);
            if (failed) failed_points_sampling.push_back(x_sample);
            int landed_rank = failed?0:part.get_rank(utils::tag[partition::global](lb));
            auto& indices = output.mine.idx[i];
            auto& coeffs  = output.mine.coeff[i];
            indices = 0;
            coeffs  = 0.0;
            
            if (!failed)
            {
                if (group.rank() == landed_rank)
                {
                    const auto block_bbx = grid.get_bounding_box(lb);
                    const auto dx = grid.get_dx(lb);
                    index_t landed_cell(0, 0, 0, part.to_local(lb));
                    ctrs::array<int, 3>     deltai      = 0;
                    ctrs::array<coeff_t, 3> reduced_idx = 0;
                    for (int d = 0; d < dim; ++d)
                    {
                        reduced_idx[d]  = (x_sample[d] - block_bbx.min(d)) / dx[d];
                        landed_cell[d]  = floor(reduced_idx[d]);
                        reduced_idx[d] -= 0.5;
                        deltai[d]       = utils::sign(reduced_idx[d] - landed_cell[d]);
                    }
                    if (!strategy.try_make_cloud(indices, coeffs, grid, landed_cell, x_sample, reduced_idx, deltai, exclusion_crit))
                    {
                        failed_points_interp.push_back(x_sample);
                    }
                }
                else
                {
                    //Compute the scatter indices here
                    off_infos.push_back(lfinfo{landed_rank, int(lb.value), x_sample, i});
                }
            }
        }
        
        std::sort(off_infos.begin(), off_infos.end(), [](const auto& a, const auto& b){ return a.rank < b.rank; });
        
        for (const auto& l: off_infos)
        {
            for (const auto& xx: l.pt) point_msg.send_buffers[l.rank].push_back(xx);
            output.scatter_offsets.push_back(l.gidx);
        }        

        constexpr bool supports_par_interp = requires { group.template send_all<int>(); };
        
        if constexpr (supports_par_interp)
        {
            point_msg.assure_recv_buf_size(group);
            point_msg.send_all(group);
            
            std::vector<pnt_t> other_points;
            for (int ii = 0; ii < group.size(); ++ii)
            {
                for (std::size_t jj = 0; jj < point_msg.recv_buffers[ii].size(); jj += pnt_t::size())
                {
                    pnt_t newpt;
                    for (int kk = 0; kk < pnt_t::size(); ++kk)
                    {
                        newpt[kk] = point_msg.recv_buffers[ii][jj + kk];
                    }
                    other_points.push_back(newpt);
                    output.send_sizes[ii]++;
                }
            }
            
            output.others.idx.resize(other_points.size());
            output.others.coeff.resize(other_points.size());
            
            for (std::size_t i = 0; i < other_points.size(); ++i)
            {
                bool failed = false;
                const auto& x_sample = other_points[i];
                const auto lb = get_glob_block(x_sample, failed);
                if (failed) failed_points_sampling.push_back(x_sample);
                int landed_rank = failed?0:part.get_rank(utils::tag[partition::global](lb));
                auto& indices = output.others.idx[i];
                auto& coeffs  = output.others.coeff[i];
                indices = 0;
                coeffs  = 0.0;
                
                if (!failed)
                {
                    if (group.rank() == landed_rank)
                    {
                        const auto block_bbx = grid.get_bounding_box(lb);
                        const auto dx = grid.get_dx(lb);
                        index_t landed_cell(0, 0, 0, part.to_local(lb));
                        ctrs::array<int, 3>     deltai      = 0;
                        ctrs::array<coeff_t, 3> reduced_idx = 0;
                        for (int d = 0; d < dim; ++d)
                        {
                            reduced_idx[d]  = (x_sample[d] - block_bbx.min(d)) / dx[d];
                            landed_cell[d]  = floor(reduced_idx[d]);
                            reduced_idx[d] -= 0.5;
                            deltai[d]       = utils::sign(reduced_idx[d] - landed_cell[d]);
                        }
                        if (!strategy.try_make_cloud(indices, coeffs, grid, landed_cell, x_sample, reduced_idx, deltai, exclusion_crit))
                        {
                            failed_points_interp.push_back(x_sample);
                        }
                    }
                    else
                    {
                        // Never get here!
                        throw except::sp_exception("A thread was sent a point for sampling that it does not contain, this is impossible!");
                    }
                }
            }
            
            constexpr int num_vars = arr_t::alias_type::size();
            for (int ii = 0; ii < group.size(); ++ii)
            {
                output.result_message.send_buffers[ii].resize(num_vars*output.send_sizes[ii]);
            }
            
            output.result_message.assure_recv_buf_size(group);
            
            //Better to fail here than later
            output.result_message.send_all(group);
        }
        
        if (failed_points_sampling.size() > 0)
        {
            throw except::points_exception("Sampling operation failed for " + std::to_string(failed_points_sampling.size()) + " points (cannot find block for sampling)", failed_points_sampling);
        }
        
        if (failed_points_interp.size() > 0)
        {
            throw except::points_exception("Sampling operation failed for " + std::to_string(failed_points_interp.size()) + " points (cannot find interp. cloud)", failed_points_interp);
        }
        
        
        output.transfer();
        
        return output;
    }
    
    template <typename arr_t, ctrs::basic_array container_t, typename strategy_t>
    static auto create_interpolation(const arr_t& arr, const container_t& points, const strategy_t& strategy)
    {
        return create_interpolation(arr, points, strategy, [](const typename arr_t::index_type&){return false;});
    }
    
    template <typename container_t, typename arr_t, typename operator_t>
    static void sample_array(container_t& output, const arr_t& arr, operator_t& oper)
    {
        output.resize(oper.mine.size());
        const auto idx_img   = utils::make_vec_image(oper.mine.idx.data(arr.device()));
        const auto coeff_img = utils::make_vec_image(oper.mine.coeff.data(arr.device()));
        
        const auto& group = arr.get_grid().group();
        
        constexpr bool supports_par_interp = requires { group.template send_all<int>(); };
        
        constexpr static bool is_shared_vec = requires { output.transfer(); };
        auto out_img = [&]()
        {
            if constexpr (is_shared_vec) return utils::make_vec_image(output.data(arr.device()));
            else                        return utils::make_vec_image(output);
        }();
        
        const auto arr_img   = arr.image();
        
        using real_t = typename arr_t::fundamental_type;
        auto kern = [=] _sp_hybrid (const std::size_t& i) mutable
        {
            auto& result = out_img[i];
            result = real_t(0.0);
            const auto& coeffs = coeff_img[i];
            const auto& idxs   = idx_img  [i];
            for (int j = 0; j < idxs.size(); ++j)
            {
                auto data = arr_img.get_elem(idxs[j]);
                result = result + coeffs[j]*data;
            }
        };
        auto rg = dispatch::ranges::from_array(out_img, arr.device());
        dispatch::execute(rg, kern);
        
        if constexpr (supports_par_interp)
        {
        
            using alias_type = typename arr_t::alias_type;
            std::size_t offset = 0;
            for (int p = 0; p < group.size(); ++p)
            {
                if (p != group.rank())
                {
                    std::size_t send_size = oper.send_sizes[p];
                    
                    const auto o_idx_img   = utils::make_vec_image(oper.others.idx.data(arr.device()));
                    const auto o_coeff_img = utils::make_vec_image(oper.others.coeff.data(arr.device()));
                    auto dest_img = utils::make_vec_image(oper.result_message.send_buffers[p].data(arr.device()));
                    
                    auto comm_sampl_kern = [=] _sp_hybrid (const std::size_t& i) mutable
                    {
                        alias_type sampled = real_t(0.0);
                        const auto& coeffs = o_coeff_img[offset + i];
                        const auto& idxs   = o_idx_img[offset + i];
                        for (int j = 0; j < idxs.size(); ++j)
                        {
                            auto data = arr_img.get_elem(idxs[j]);
                            sampled = sampled + coeffs[j]*data;
                        }
                        std::size_t base = i*alias_type::size();
                        for (int v = 0; v < alias_type::size(); ++v)
                        {
                            dest_img[base+v] = sampled[v];
                        }
                    };
                    
                    auto smallrg  = dispatch::ranges::make_range(0UL, send_size);
                    
                    dispatch::execute(smallrg, comm_sampl_kern, arr.device());
                    
                    offset += send_size;
                }
            }
            
            oper.result_message.send_all(group);
            
            std::size_t scatter_offst = 0;
            for (int p = 0; p < group.size(); ++p)
            {
                if (p != group.rank())
                {
                    auto buf_img     = utils::make_vec_image(oper.result_message.recv_buffers[p].data(arr.device()));
                    auto scatter_img = utils::make_vec_image(oper.scatter_offsets.data(arr.device()));
                    auto scatter_kern = [=] _sp_hybrid (const std::size_t& i) mutable
                    {
                        alias_type result;
                        std::size_t base = i*alias_type::size();
                        for (int v = 0; v < alias_type::size(); ++v)
                        {
                            result[v] = buf_img[base + v];
                        }
                        std::size_t scatter_idx = scatter_img[scatter_offst + i];
                        out_img[scatter_idx] = result;
                    };
                    std::size_t num_elems = oper.result_message.recv_buffers[p].size() / alias_type::size();
                    scatter_offst += num_elems;
                    auto smallrg  = dispatch::ranges::make_range(0UL, num_elems);
                    dispatch::execute(smallrg, scatter_kern, arr.device());
                }
            }
        }
        
        if constexpr (is_shared_vec && device::is_gpu<typename arr_t::device_type>) output.itransfer();
    }
    
    //todo: update this to include an omni kernel
    template <typename arr_t, typename operator_t>
    [[nodiscard]] static auto sample_array(const arr_t& arr, operator_t& oper)
    {
        using alias_type = typename arr_t::alias_type;
        using cpu_out_t  = std::vector<alias_type>;
        using gpu_out_t  = device::device_vector<alias_type>;
        constexpr static bool is_gpu = device::is_gpu<typename arr_t::device_type>;
        using output_type = std::conditional<is_gpu, gpu_out_t, cpu_out_t>::type;
        output_type output;
        sample_array(output, arr, oper);
        return output;
    }
}
