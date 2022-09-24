#pragma once

#include "core/grid.h"
#include "fetch/fetch.h"

namespace spade::pde_algs
{    
    template <typename flux_func_t> concept is_flux_functor =
    requires(const flux_func_t& f)
    {
        //TODO
        f;
    };
    
    template <typename arr_t, typename flux_eval_t> concept has_flux_compatibility =
    requires(const arr_t& t0, const flux_eval_t& t1)
    {
        //TODO
        t0;
    };
    
    enum block_flux_policy
    {
        block_flux_all,
        block_flux_interior,
        block_flux_boundary
    };
    
    
    
    template <
        grid::multiblock_array array1_t,
        grid::multiblock_array array2_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<array1_t, grid::cell_centered> &&
        (dims::rank_eval<typename array1_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename array1_t::array_major_dim_t>::value==0)
    static void flux_div(
        const array1_t& prims,
        array2_t& rhs,
        const block_flux_policy& policy,
        const ctrs::array<bool, array1_t::grid_type::dim()>& domain_boundary_flux,
        const flux_funcs_t&... flux_funcs)
    {
        typedef typename array1_t::value_type real_type;
        const grid::multiblock_grid auto& ar_grid = prims.get_grid();
        int i3d = ((ar_grid.dim()==3)?1:0);
        auto grid_range = range(-1,ar_grid.get_num_cells(0))*range(-1,ar_grid.get_num_cells(1))*range(-i3d,ar_grid.get_num_cells(2))*range(0,ar_grid.get_num_local_blocks());
        for (auto idx: grid_range*range(0, ar_grid.dim()))
        {
            int idir = idx[4];
            const real_type dx = ar_grid.get_dx(idir);
            grid::cell_idx_t il(idx[0], idx[1], idx[2], idx[3]);
            grid::cell_idx_t ir(idx[0], idx[1], idx[2], idx[3]);
            ir[idir] += 1;
            grid::face_idx_t iface = grid::cell_to_face(il, idir, 1);
            const ctrs::array<real_type,3> xyz_comp_l = ar_grid.get_comp_coords(il);
            const ctrs::array<real_type,3> xyz_comp_r = ar_grid.get_comp_coords(ir);
            const real_type jac_l = coords::calc_jacobian(ar_grid.coord_sys(), xyz_comp_l, il);
            const real_type jac_r = coords::calc_jacobian(ar_grid.coord_sys(), xyz_comp_r, ir);
            utils::foreach_param([&](const auto& flux_func) -> void
            {
                using flux_func_t = decltype(flux_func);
                typename utils::remove_all<flux_func_t>::type::input_type flux_data;
                fetch::get_flux_data(ar_grid, prims, iface, flux_data);
                auto flux = flux_func.calc_flux(flux_data);
                for (int n = 0; n < flux.size(); ++n)
                {
                    rhs(n, il[0], il[1], il[2], idx[3])-=jac_l*flux[n]/(dx);
                    rhs(n, ir[0], ir[1], ir[2], idx[3])+=jac_r*flux[n]/(dx);
                }
            }, flux_funcs...);
        }
    }
    
    template <
        grid::multiblock_array array1_t,
        grid::multiblock_array array2_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<array1_t, grid::cell_centered> &&
        (dims::rank_eval<typename array1_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename array1_t::array_major_dim_t>::value==0)
    static void flux_div(
        const array1_t& prims,
        array2_t& rhs,
        const flux_funcs_t&... flux_funcs)
        {
            const block_flux_policy policy = block_flux_all;
            const ctrs::array<bool, array1_t::grid_type::dim()> domain_boundary_flux(true);
            flux_div(prims, rhs, policy, domain_boundary_flux, flux_funcs...);
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
