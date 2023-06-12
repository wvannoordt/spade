#pragma once

#include "grid/grid.h"
#include "omni/omni.h"

// GET RID OF THIS IMPLEMENTATION, ITS EXACTLY THE SAME AS FLUX_DIV!!!!!!!!

namespace spade::pde_algs
{
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void boundary_flux(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const bound_box_t<bool,sol_arr_t::grid_type::dim()>& apply_to_boundary,
        const flux_funcs_t&... flux_funcs)
    {

        using omni_union_type = omni::combine_omni_stencils<flux_funcs_t...>;
        using alias_type = sol_arr_t::alias_type;
        const auto& grid = prims.get_grid();
        for (auto lb: range(0, grid.get_num_local_blocks()))
        {
            const auto& lb_glob = grid.get_partition().get_global_block(lb);
            const auto& idomain = grid.is_domain_boundary(lb_glob);
            for (int dir = 0; dir < 2*grid.dim(); ++dir)
            {
                int rdir = dir/2;
                int mdir = dir%2;
                const bool apply_to_cur_boundary = idomain(rdir, mdir) && apply_to_boundary(rdir, mdir);
                if (apply_to_cur_boundary)
                {
                    bound_box_t<int,3> flux_bounds;
                    flux_bounds.min(0) = -1;
                    flux_bounds.min(1) = -1;
                    flux_bounds.min(2) = -1;
                    if constexpr (grid.dim()==2) flux_bounds.min(2) = 0;

                    flux_bounds.max(0) = grid.get_num_cells(0);
                    flux_bounds.max(1) = grid.get_num_cells(1);
                    flux_bounds.max(2) = grid.get_num_cells(2);
                    const auto dx = grid.get_dx(rdir);
                    if (mdir == 0)
                    {
                        flux_bounds.max(rdir) = 0;
                    }
                    else
                    {
                        flux_bounds.min(rdir) = grid.get_num_cells(rdir)-1;
                    }
                    auto g0 = range(flux_bounds.min(0),flux_bounds.max(0)); //x
                    auto g1 = range(flux_bounds.min(1),flux_bounds.max(1)); //y
                    auto g2 = range(flux_bounds.min(2),flux_bounds.max(2)); //z
                    auto boundary_range = g0*g1*g2;
                    for (auto idx: boundary_range)
                    {
                        const auto dx = grid.get_dx(rdir);
                        grid::cell_idx_t il(idx[0_c], idx[1_c], idx[2_c], lb);
                        grid::cell_idx_t ir(idx[0_c], idx[1_c], idx[2_c], lb);
                        ir[rdir] += 1;
                        grid::face_idx_t iface = grid::cell_to_face(il, rdir, 1);
                        const auto xyz_comp_l = grid.get_comp_coords(il);
                        const auto xyz_comp_r = grid.get_comp_coords(ir);
                        const auto jac_l = coords::calc_jacobian(grid.coord_sys(), xyz_comp_l, il);
                        const auto jac_r = coords::calc_jacobian(grid.coord_sys(), xyz_comp_r, ir);
                        
                        using omni_type = omni_union_type;
                        using data_type = omni::stencil_data_t<omni_type, sol_arr_t>;
                        data_type data;
                        omni::retrieve(prims, iface, data);
                        using flux_out_t = rhs_arr_t::alias_type;
                        flux_out_t accum(0.0);
                        utils::foreach_param([&](const auto& flux_func)
                        {   
                            using flux_func_t = decltype(flux_func);
                            using omni_type   = utils::remove_all<flux_func_t>::type::omni_type;

                            auto flux_data = omni::interpret_stencil<omni_type>(data);
                            accum += flux_func(flux_data);
                            
                        }, flux_funcs...);

                        //todo: fix this garbage
                        if constexpr (ctrs::basic_array<alias_type>)
                        {
                            for (int n = 0; n < accum.size(); ++n)
                            {
                                rhs(n, il) -= jac_l*accum[n]/(dx);
                                rhs(n, ir) += jac_r*accum[n]/(dx);
                            }
                        }
                        else
                        {
                            rhs(il) -= jac_l*accum/(dx);
                            rhs(ir) += jac_r*accum/(dx);
                        }
                    }
                }
            }
        }
    }
}