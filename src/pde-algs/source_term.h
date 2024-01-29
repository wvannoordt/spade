#pragma once

#include <concepts>

#include "grid/grid.h"
#include "omni/omni.h"

namespace spade::pde_algs
{
    namespace detail
    {
        template <grid::multiblock_array sol_arr_t, typename source_term_t>
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
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void source_term(const sol_arr_t& q, rhs_arr_t& rhs, const source_term_t& source_term_func)
    {
        typedef typename sol_arr_t::value_type real_type;
        const grid::multiblock_grid auto& ar_grid = rhs.get_grid();
        const grid::array_centering ctr = sol_arr_t::centering_type();
        const auto kernel = omni::to_omni<ctr>(source_term_func, q);
        const auto qv = q.image();
        auto rv       = rhs.image();

        const auto grid_img = ar_grid.image(q.device());
        auto grid_range = dispatch::support_of(q, grid::exclude_exchanges);
        auto loop = [=] _sp_hybrid (const grid::cell_idx_t& icell) mutable
        {
            const auto xc       = grid_img.get_comp_coords(icell);
            const real_type jac = coords::calc_jacobian(grid_img.get_coord_sys(), xc, icell);
            using kernel_t      = decltype(kernel); 
            using omni_type     = utils::remove_all<kernel_t>::type::omni_type;
            using input_type    = omni::stencil_data_t<omni_type, sol_arr_t>;

            input_type source_data;
            omni::retrieve(grid_img, qv, icell, source_data);
            auto source_term = kernel(source_data);
            rv.incr_elem(icell, source_term/jac);
        };
        
        dispatch::execute(grid_range, loop);
    }
}