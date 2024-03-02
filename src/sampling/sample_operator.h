#pragma once

#include "grid/grid_index_types.h"
#include "sampling/sample_clouds.h"
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
        
        const auto& grid = arr.get_grid();
        constexpr int dim = grid.dim();
        
        // We will use a BVH to check which block each point is in.
        // Might want to consider pre-computing this somewhere
        geom::bvh_t<dim, typename pnt_t::value_type> block_bvh;
        geom::bvh_params_t bvhparams{4, 1000};
        
        // Perform in computational coordinates, note that we strictly
        // need a 2D BVH if the grid is 2D
        const auto bbx = grid.get_bounds();
        using bnd_t = bound_box_t<typename pnt_t::value_type, dim>;
        bnd_t bnd;
        for (int d = 0; d < dim; ++d)
        {
            bnd.min(d) = bbx.min(d);
            bnd.max(d) = bbx.max(d);
        }
        bnd = bnd.inflate(1.1);
        
        const auto el_check = [&](const std::size_t& lb_glob, const auto& bnd_in)
        {
            const auto lb = utils::tag[partition::global](lb_glob);
            const auto block_box = grid.get_bounding_box(lb).inflate(1.1);
            return block_box.intersects(bnd_in);
        };
        
        block_bvh.calculate(bnd, grid.get_num_global_blocks(), el_check, bvhparams);

        //We will now compute the indices of each point in the interpolation cloud using the BVH
        constexpr int num_coeffs = strategy_t::stencil_size();
        using output_t = interp_oper_t<grid_t, coeff_t, arr_t::centering_type(), num_coeffs>;
        output_t output;
        
        output.idx.resize(points.size());
        output.coeff.resize(points.size());
        
        std::vector<pnt_t> failed_points;
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            const auto& x_sample = points[i];
            auto lb = utils::tag[partition::global](-1);
            const auto eval = [&](const std::size_t& lb_cand)
            {
                auto lb_tagged = utils::tag[partition::global](lb_cand);
                const auto block_box = grid.get_bounding_box(lb_tagged);
                if (block_box.contains(x_sample)) lb.value = lb_cand;
            };
            
            using vec_t = typename decltype(block_bvh)::pnt_t;
            block_bvh.check_elements(eval, ctrs::to_array(x_sample));
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
                    std::stringstream ss;
                    ss << std::setprecision(20);
                    ss << "Cannot compute find suitable block for point:\n" << x_sample;
                    ss << "\nBounds:\n" << bbx;
                    ss << "\nNum. checks: " << num_checked;
                    std::vector<pnt_t> points{x_sample};
                    throw except::points_exception<typename pnt_t::value_type>("Error attempting to build sampling operator:\n" + ss.str(), points);
                }
            }
            
            const auto block_bbx = grid.get_bounding_box(lb);
            const auto dx = grid.get_dx(lb);
            
            const auto& partition = grid.get_partition();
            using index_t = output_t::index_type;
            index_t landed_cell(0, 0, 0, partition.to_local(lb));
            if (landed_cell.lb() == partition.no_value)
            {
                throw except::sp_exception("cannot currently perform cross-process interpolation");
            }
            
            ctrs::array<int, 3>     deltai      = 0;
            ctrs::array<coeff_t, 3> reduced_idx = 0;
            for (int d = 0; d < dim; ++d)
            {
                reduced_idx[d]  = (x_sample[d] - block_bbx.min(d)) / dx[d];
                landed_cell[d]  = floor(reduced_idx[d]);
                reduced_idx[d] -= 0.5;
                deltai[d] = utils::sign(reduced_idx[d] - landed_cell[d]);
            }
            
            auto& indices = output.idx[i];
            auto& coeffs  = output.coeff[i];
            
            indices = landed_cell;
            coeffs  = 0.0;
            
            if (!strategy.try_make_cloud(indices, coeffs, grid, landed_cell, x_sample, reduced_idx, deltai, exclusion_crit))
            {
                failed_points.push_back(x_sample);
            }
        }
        
        if (failed_points.size() > 0)
        {
            throw except::points_exception("Sampling operation failed for " + std::to_string(failed_points.size()) + " points", failed_points);
        }
        
        output.idx.transfer();
        output.coeff.transfer();
        return output;
    }
    
    template <typename arr_t, ctrs::basic_array container_t, typename strategy_t>
    static auto create_interpolation(const arr_t& arr, const container_t& points, const strategy_t& strategy)
    {
        return create_interpolation(arr, points, strategy, [](const typename arr_t::index_type&){return false;});
    }
    
    template <typename container_t, typename arr_t, typename operator_t>
    static void sample_array(container_t& output, const arr_t& arr, const operator_t& oper)
    {
        output.resize(oper.size());
        const auto idx_img   = utils::make_vec_image(oper.idx.data(  arr.device()));
        const auto coeff_img = utils::make_vec_image(oper.coeff.data(arr.device()));
        
        constexpr static bool is_shared_vec = requires { output.transfer(); };
        auto out_img   = [&]()
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
        
        if constexpr (is_shared_vec && device::is_gpu<typename arr_t::device_type>) output.itransfer();
    }
    
    //todo: update this to include an omni kernel
    template <typename arr_t, typename operator_t>
    [[nodiscard]] static auto sample_array(const arr_t& arr, const operator_t& oper)
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
