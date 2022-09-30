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
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered> &&
        (dims::rank_eval<typename sol_arr_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename sol_arr_t::array_major_dim_t>::value==0)
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const block_flux_policy& policy,
        const bound_box_t<bool,sol_arr_t::grid_type::dim()>& domain_boundary_flux,
        const flux_funcs_t&... flux_funcs)
    {
        typedef typename sol_arr_t::value_type real_type;
        const grid::multiblock_grid auto& ar_grid = prims.get_grid();
        
        auto block_range = range(0,ar_grid.get_num_local_blocks());
        for (auto lb_idir: block_range*range(0, ar_grid.dim()))
        {
            int lb_loc  = lb_idir[0];
            int lb_glob = ar_grid.get_partition().get_global_block(lb_loc);
            int idir = lb_idir[1];
            bound_box_t<int,3> flux_bounds;
            for (auto d:range(0,3))
            {
                flux_bounds.min(d) = (ar_grid.dim()==2 && d==2)?0:-1;
                flux_bounds.max(d) = ar_grid.get_num_cells(d);
            }
            const auto& iboundary = ar_grid.is_domain_boundary(lb_glob);
            for (auto d:range(0,ar_grid.dim()))
            {
                if (iboundary.min(d) && !domain_boundary_flux.min(d)) flux_bounds.min(d) = 0;
                if (iboundary.max(d) && !domain_boundary_flux.max(d)) flux_bounds.max(d) = ar_grid.get_num_cells(d) - 1;
            }
            auto g0 = range(flux_bounds.min(0),flux_bounds.max(0));
            auto g1 = range(flux_bounds.min(1),flux_bounds.max(1));
            auto g2 = range(flux_bounds.min(2),flux_bounds.max(2));
            
            auto grid_range = g0*g1*g2;
            for (auto idx: grid_range)
            {
                const real_type dx = ar_grid.get_dx(idir);
                grid::cell_idx_t il(idx[0], idx[1], idx[2], lb_loc);
                grid::cell_idx_t ir(idx[0], idx[1], idx[2], lb_loc);
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
                    fetch::get_face_data(ar_grid, prims, iface, flux_data);
                    auto flux = flux_func.calc_flux(flux_data);
                    for (int n = 0; n < flux.size(); ++n)
                    {
                        rhs(n, il[0], il[1], il[2], il[3]) -= jac_l*flux[n]/(dx);
                        rhs(n, ir[0], ir[1], ir[2], ir[3]) += jac_r*flux[n]/(dx);
                    }
                }, flux_funcs...);
            }
        }
    }
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered> &&
        (dims::rank_eval<typename sol_arr_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename sol_arr_t::array_major_dim_t>::value==0)
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_funcs_t&... flux_funcs)
    {
        const std::size_t dim = sol_arr_t::grid_type::dim();
        const block_flux_policy policy = block_flux_all;
        bound_box_t<bool,dim> domain_boundary_flux;
        for (auto i: range(0,dim))
        {
            domain_boundary_flux.max(i) = true;
            domain_boundary_flux.min(i) = true;
        }
        flux_div(prims, rhs, policy, domain_boundary_flux, flux_funcs...);
    }

    template <typename T> concept variate_source_term = fetch::fetch_type<typename T::input_type, fetch::fetch_cell> && requires (T t)
    {
        typename T::input_type;
        typename T::output_type;
    };
    
    namespace detail
    {
        template <grid::multiblock_array sol_arr_t, typename source_term_t>
        requires (variate_source_term<source_term_t>)
        auto eval_source_term(const grid::cell_idx_t& ijk, const sol_arr_t& q, const source_term_t& source_term_func)
        {
            using out_t = typename source_term_t::output_type;
            using in_t  = typename source_term_t::input_type;
            in_t input;
            fetch::get_cell_data(q.get_grid(), q, ijk, input);
            out_t output = source_term_func(input);
            return output;
        }
        
        template <grid::multiblock_array sol_arr_t, typename source_term_t>
        requires (!variate_source_term<source_term_t>)
        auto eval_source_term(const grid::cell_idx_t& ijk, const sol_arr_t& q, const source_term_t& source_term_func)
        {
            return source_term_func();
        }
    }
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename source_term_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered> &&
        (dims::rank_eval<typename sol_arr_t::array_minor_dim_t>::value==1) &&
        (dims::rank_eval<typename sol_arr_t::array_major_dim_t>::value==0)
    static void source_term(const sol_arr_t& q, rhs_arr_t& rhs, const source_term_t& source_term_func)
    {
        typedef typename sol_arr_t::value_type real_type;
        const grid::multiblock_grid auto& ar_grid = rhs.get_grid();
        auto grid_range = ar_grid.get_range(sol_arr_t::centering_type(), grid::exclude_exchanges);
        for (auto idx: grid_range)
        {
            ctrs::array<grid::cell_t<int>, 4> i(idx[0], idx[1], idx[2], idx[3]);
            const ctrs::array<real_type,3> xc = ar_grid.get_comp_coords(i);
            const real_type jac = coords::calc_jacobian(ar_grid.coord_sys(), xc, i);
            
            const ctrs::basic_array auto source_term = detail::eval_source_term(i, q, source_term_func);
            for (auto n: range(0, source_term.size()))
            {
                rhs(n, i[0], i[1], i[2], i[3])+=source_term[n]/jac;
            }
        }
    }
}
