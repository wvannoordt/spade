#pragma once

#include "algs/alg_traits.h"
#include "pde-algs/flux-div/tags.h"
#include "pde-algs/flux-div/flux_div_basic.h"
#include "pde-algs/flux-div/flux_div_ldbal.h"
#include "pde-algs/flux-div/flux_div_bfoct.h"
#include "pde-algs/flux-div/flux_div_fldbc.h"

namespace spade::pde_algs
{
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t,
        typename alg_traits_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func,
        const alg_traits_t& traits)
    {
        using namespace sym::literals;
        const auto flux_algorithm = algs::get_trait(traits, "flux_div"_sym);
        using fdiv_tag_t = typename utils::remove_all<decltype(flux_algorithm)>::type;
        if constexpr (std::same_as<fdiv_tag_t, tbasic_t>       ) flux_div_basic(prims, rhs, flux_func, traits);
        if constexpr (std::same_as<fdiv_tag_t, tbfoct_t>       ) flux_div_bfoct(prims, rhs, flux_func, traits);
        if constexpr (std::same_as<fdiv_tag_t, tldbal_t<true>> ) flux_div_ldbal(prims, rhs, flux_func, tldbal_t<true>(),  traits);
        if constexpr (std::same_as<fdiv_tag_t, tldbal_t<false>>) flux_div_ldbal(prims, rhs, flux_func, tldbal_t<false>(), traits);
        if constexpr (std::same_as<fdiv_tag_t, tfldbc_t<true>> ) flux_div_fldbc(prims, rhs, flux_func, tfldbc_t<true>(),  traits);
        if constexpr (std::same_as<fdiv_tag_t, tfldbc_t<false>>) flux_div_fldbc(prims, rhs, flux_func, tfldbc_t<false>(), traits);
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
        const auto default_traits = algs::make_traits(basic);
        flux_div(prims, rhs, flux_func, default_traits);
    }
}


