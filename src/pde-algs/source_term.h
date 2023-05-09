#pragma once

#include <concepts>

#include "core/grid.h"
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

        auto grid_range = ar_grid.get_range(sol_arr_t::centering_type(), grid::exclude_exchanges);
        for (auto idx: grid_range)
        {
            grid::cell_idx_t i;
            i.i()  = idx[0];
            i.j()  = idx[1];
            i.k()  = idx[2];
            i.lb() = idx[3];
            const auto xc = ar_grid.get_comp_coords(i);
            const real_type jac = coords::calc_jacobian(ar_grid.coord_sys(), xc, i);
            using kernel_t    = decltype(kernel);
            using omni_type   = utils::remove_all<kernel_t>::type::omni_type;
            using input_type  = omni::stencil_data_t<omni_type, sol_arr_t>;

            input_type source_data;
            omni::retrieve(q, i, source_data);
            auto source_term = kernel(source_data);
            for (auto n: range(0, source_term.size()))
            {
                rhs(n, i[0], i[1], i[2], i[3])+=source_term[n]/jac;
            }
        }
    }
}