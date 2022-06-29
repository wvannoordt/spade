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
        return 0.0;
    }
}