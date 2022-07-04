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
        const int idir = boundary_id/2;
        //1 = positive face, 0 = negative
        const int pm   = boundary_id%2;
        spade::bound_box_t<int,3> bnd;
        for (auto d: range(0,3))
        {
            bnd.min(d) = 0;
            bnd.max(d) = grid.get_num_cells(d);
        }
        bnd.min(idir) = pm*grid.get_num_cells(idir);
        bnd.max(idir) = pm*grid.get_num_cells(idir);
        auto flx_rg = range(bnd.min(0), bnd.max(0))*range(bnd.min(1), bnd.max(1))*range(bnd.min(2), bnd.max(2));
        typedef typename flux_func_t::output_type out_t;
        out_t ini = 0.0;
        spade::reduce_ops::reduce_sum<out_t> sum_op;
        sum_op.init(ini);
        for (auto lb: range(0, grid.get_num_local_blocks()))
        {
            if (grid.is_domain_boundary(lb, boundary_id))
            {
                for (auto i: flx_rg)
                {
                    spade::grid::cell_idx_t i_c (i[0], i[1], i[2], lb);
                    spade::grid::face_idx_t i_f = spade::grid::cell_to_face(i_c, idir, pm);
                }
            }
        }
        return 0.0;
    }
}