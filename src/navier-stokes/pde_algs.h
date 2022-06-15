#pragma once

#include "core/grid.h"

#include "navier-stokes/convective.h"

namespace cvdf::pde_algs
{
    template <
        grid::multiblock_array array1_t,
        grid::multiblock_array array2_t,
        typename flux_func_t>
    requires
        grid::has_centering_type<array1_t, grid::cell_centered> &&
        (dims::rank_eval<typename array1_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename array1_t::array_major_dim_t>::value==0)
    static void flux_div(const array1_t& prims, array2_t& rhs, const flux_func_t& flux_func)
    {
        typedef typename array1_t::value_type real_type;
        const grid::multiblock_grid auto& ar_grid = prims.get_grid();
        int i3d = ((ar_grid.dim()==3)?1:0);
        auto grid_range = range(-1,ar_grid.get_num_cells(0))*range(-1,ar_grid.get_num_cells(1))*range(-i3d,ar_grid.get_num_cells(2))*range(0,ar_grid.get_num_local_blocks());
        for (auto idx: grid_range*range(0, ar_grid.dim()))
        {
            int idir = idx[4];
            ctrs::array<grid::face_t<int>, 5> iface(idir, idx[0], idx[1], idx[2], idx[3]);
            ctrs::array<grid::cell_t<int>, 4> il(idx[0], idx[1], idx[2], idx[3]);
            ctrs::array<grid::cell_t<int>, 4> ir(idx[0], idx[1], idx[2], idx[3]);
            ir[idir] += 1;
            const ctrs::array<real_type,3> xyz_comp_l = ar_grid.get_comp_coords(il[0], il[1], il[2], il[3]);
            const ctrs::array<real_type,3> xyz_comp_r = ar_grid.get_comp_coords(ir[0], ir[1], ir[2], ir[3]);
            typename flux_func_t::input_type flux_data;
            flux_input::get_flux_data(ar_grid, prims, iface, flux_data);
            const real_type jac_l = coords::calc_jacobian(ar_grid.coord_sys(), xyz_comp_l, il);
            const real_type jac_r = coords::calc_jacobian(ar_grid.coord_sys(), xyz_comp_r, ir);
            fluid_state::flux_t<real_type> flux = flux_func.calc_flux(flux_data);
            const real_type dx = ar_grid.get_dx(idir);
            for (int n = 0; n < flux.size(); ++n)
            {
                rhs(n, il[0], il[1], il[2], idx[3])-=jac_l*flux[n]/(dx);
                rhs(n, ir[0], ir[1], ir[2], idx[3])+=jac_r*flux[n]/(dx);
            }
        }
    }
    
    template <
        grid::multiblock_array array_t,
        typename source_term_t>
    requires
        grid::has_centering_type<array_t, grid::cell_centered> &&
        (dims::rank_eval<typename array_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename array_t::array_major_dim_t>::value==0)
    static void source_term(array_t& rhs, const source_term_t& source_term_func)
    {
        typedef typename array_t::value_type real_type;
        const grid::multiblock_grid auto& ar_grid = rhs.get_grid();
        auto grid_range = ar_grid.get_range(array_t::centering_type(), grid::exclude_exchanges);
        for (auto idx: grid_range)
        {
            ctrs::array<grid::cell_t<int>, 4> i(idx[0], idx[1], idx[2], idx[3]);
            const ctrs::array<real_type,3> xc = ar_grid.get_comp_coords(i);
            const real_type jac = coords::calc_jacobian(ar_grid.coord_sys(), xc, i);
            //TODO: generalize this
            const ctrs::basic_array auto source_term = source_term_func();
            for (auto n: range(0, source_term.size()))
            {
                rhs(n, i[0], i[1], i[2], i[3])+=source_term[n]/jac;
            }
        }
    }
}
