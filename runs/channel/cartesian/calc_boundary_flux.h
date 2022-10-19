#include "spade.h"

#pragma once

namespace proto
{
    template <spade::grid::multiblock_array array_t, spade::pde_algs::is_flux_functor flux_func_t>
    requires spade::pde_algs::has_flux_compatibility<array_t, flux_func_t> &&
        spade::grid::has_centering_type<array_t, spade::grid::cell_centered> &&
        (spade::dims::rank_eval<typename array_t::array_minor_dim_t>::value==1) &&
        (spade::dims::rank_eval<typename array_t::array_major_dim_t>::value==0)
    static auto get_domain_boundary_flux(const array_t& prim, const flux_func_t& func, const int& boundary_id)
    {
        const auto& grid = prim.get_grid();
        typedef typename array_t::grid_type::coord_type coord_type;
        coord_type dV = 1.0;
        for (auto d: range(0,grid.dim())) dV *= grid.get_dx(d);
        const int idir = boundary_id/2;
        const coord_type dA = dV/grid.get_dx(idir);
        //1 = positive face, 0 = negative
        const int pm   = boundary_id%2;
        spade::bound_box_t<int,3> bnd;
        for (auto d: range(0,3))
        {
            bnd.min(d) = 0;
            bnd.max(d) = grid.get_num_cells(d);
        }
        bnd.min(idir) = pm*(grid.get_num_cells(idir)-1);
        bnd.max(idir) = pm*(grid.get_num_cells(idir)-1);
        auto flx_rg = range(bnd.min(0), bnd.max(0))*range(bnd.min(1), bnd.max(1))*range(bnd.min(2), bnd.max(2));
        typedef typename flux_func_t::output_type out_t;
        out_t ini = 0.0;
        spade::reduce_ops::reduce_sum<out_t> sum_op;
        sum_op.init(ini);
        for (auto lb: range(0, grid.get_num_local_blocks()))
        {
            auto lb_glob = grid.get_partition().get_global_block(lb);
            const auto& idomain = grid.is_domain_boundary(lb_glob);
            if (idomain(boundary_id/2, boundary_id%2))
            {
                for (auto i: flx_rg)
                {
                    spade::grid::cell_idx_t i_c (i[0], i[1], i[2], lb);
                    spade::grid::face_idx_t i_f = spade::grid::cell_to_face(i_c, idir, pm);
                    auto x = grid.get_coords(i_f);
                    typename flux_func_t::input_type flux_data;
                    spade::fetch::get_face_data(grid, prim, i_f, flux_data);
                    // if (grid.group().isroot() && boundary_id == 3)
                    // {
                    //     const auto& ql = std::get<0>(flux_data.cell_data.left.elements).data;
                    //     const auto& qr = std::get<0>(flux_data.cell_data.right.elements).data;
                    //     print(i_c, i_f);
                    //     print(spade::grid::face_to_cell(i_f, 0));
                    //     print(spade::grid::face_to_cell(i_f, 1));
                    //     print("ql", ql);
                    //     print("qr", qr);
                    //     print("flx", func.calc_flux(flux_data));
                    //     grid.group().pause();
                    // }
                    auto flux = func.calc_flux(flux_data);
                    flux *= dA;
                    sum_op.reduce_elem(flux);
                }
            }
        }
        out_t loc = sum_op.value;
        out_t output = 0.0;
        for (auto i: range(0, output.size())) output[i] = grid.group().sum(loc[i]);
        return output;
    }
}
