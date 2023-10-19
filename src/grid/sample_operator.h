#pragma once

#include "grid/grid_index_types.h"
#include "core/ctrs.h"

namespace spade::grid
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
    
    template <typename arr_t, ctrs::basic_array container_t>
    static auto create_interpolation(const arr_t& arr, const container_t& points)
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
        
        // Perform in computational coordinates, node that we strictly
        // need a 2D BVH ig the grid is 2D
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
            const auto block_box = grid.get_bounding_box(lb);
            return block_box.intersects(bnd_in);
        };
        
        block_bvh.calculate(bnd, grid.get_num_global_blocks(), el_check);

        //We will now compute the indices of each point in the interpolation cloud using the BVH
        using output_t = interp_oper_t<grid_t, coeff_t, arr_t::centering_type(), 12>;
        output_t output;
        
        output.idx.resize(points.size());
        output.coeff.resize(points.size());
        
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            const auto& x_sample = points[i];
            auto lb = utils::tag[partition::global](-1);
            const auto eval = [&](const std::size_t& lb_cand)
            {
                const auto block_box = grid.get_bounding_box(utils::tag[partition::global](lb_cand));
                if (block_box.contains(x_sample)) lb.value = lb_cand;
            };
            
            using vec_t = typename decltype(block_bvh)::pnt_t;
            block_bvh.check_elements(eval, ctrs::to_array(x_sample));
            
            if (lb.value < 0)
            {
                print(lb.value, x_sample);
                print(bbx);
                print("sadyy", __FILE__, __LINE__);
                std::cin.get();
            }
            const auto block_bbx = grid.get_bounding_box(lb);
            const auto dx = grid.get_dx(lb);
            
            const auto& partition = grid.get_partition();
            using index_t = output_t::index_type;
            index_t landed_cell(0,0,0,partition.to_local(lb));
            if (landed_cell.lb() == partition.no_value)
            {
                throw except::sp_exception("cannot currently perform cross-process interpolation");
            }
            
            ctrs::array<int, 3>     deltai      = 0;
            ctrs::array<coeff_t, 3> reduced_idx = 0;
            for (int d = 0; d < dim; ++d)
            {
                reduced_idx[d] = (x_sample[d] - block_bbx.min(d)) / dx[d];
                landed_cell[d] = floor(reduced_idx[d]);
                deltai[d] = utils::sign(reduced_idx[d] - (landed_cell[d] + 0.5));
            }
            
            auto& indices = output.idx[i];
            auto& coeffs  = output.coeff[i];
            
            indices = landed_cell;
            coeffs  = 0.0;
            
            int nmx = 4;
            if (dim==3) nmx = 8;
            for (int ct = 0; ct < nmx; ++ct)
            {
                ctrs::array<int, 3> di = 0;
                di[0] = ct%2;
                di[1] = (ct/2)%2;
                di[2] = (ct/4)%2;
                
                indices[ct].i() = landed_cell.i() + di[0]*deltai[0];
                indices[ct].j() = landed_cell.j() + di[1]*deltai[1];
                indices[ct].k() = landed_cell.k() + di[2]*deltai[2];
                
                //Need to check these!
                coeff_t coeff = 1.0;
                for (int d = 0; d < dim; ++d)
                {
                    coeff_t dir_coeff0 = reduced_idx[d] - landed_cell[d];
                    coeff_t dir_coeff1 = coeff_t(1.0) - (reduced_idx[d]-landed_cell[d]);
                    coeff_t dir_coeff  = di[d]*dir_coeff0 + (coeff_t(1.0)-di[d])*dir_coeff1;
                    coeff *= dir_coeff;
                }
                coeffs[ct] = coeff;
            }
        }
        
        output.idx.transfer();
        output.coeff.transfer();
        return output;
    }
    
    template <typename container_t, typename arr_t, typename operator_t>
    static void sample_array(container_t& output, const arr_t& arr, const operator_t& oper)
    {
        output.resize(oper.size());
        const auto idx_img   = utils::make_vec_image(oper.idx.data(  arr.device()));
        const auto coeff_img = utils::make_vec_image(oper.coeff.data(arr.device()));
        auto       out_img   = utils::make_vec_image(output);
        const auto arr_img   = arr.image();
        
        using real_t = typename arr_t::fundamental_type;
        
        auto kern = _sp_lambda (const std::size_t& i) mutable
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
