#pragma once

#include "pde-algs/flux-div/tags.h"
#include "pde-algs/flux-div/flux_div_basic.h"
#include "pde-algs/flux-div/flux_div_ldbal.h"
#include "pde-algs/flux-div/flux_div_bfoct.h"

namespace spade::pde_algs
{
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t,
        typename fdiv_tag_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func,
        const fdiv_tag_t& tag)
    {
        if constexpr (std::same_as<fdiv_tag_t, tbasic_t>       ) flux_div_basic(prims, rhs, flux_func);
        if constexpr (std::same_as<fdiv_tag_t, tbfoct_t>       ) flux_div_bfoct(prims, rhs, flux_func);
        if constexpr (std::same_as<fdiv_tag_t, tldbal_t<true>> ) flux_div_ldbal(prims, rhs, flux_func, tag);
        if constexpr (std::same_as<fdiv_tag_t, tldbal_t<false>>) flux_div_ldbal(prims, rhs, flux_func, tag);
    }
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func)
    {
        flux_div(prims, rhs, flux_func, basic);
    }
}


