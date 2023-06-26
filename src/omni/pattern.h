#pragma once

#include "omni/stencil_union.h"

namespace spade::omni
{
    template <typename x_init_t, typename x_incr_t, const int remaining, const int num_points, typename infos_t, typename accum_t>
    requires (num_points > 0)
    struct pattern_impl_t
    {
        constexpr static int index = num_points - remaining;
        using xi_t = offset_t<
            x_init_t::del_i() + index*x_incr_t::del_i(),
            x_init_t::del_j() + index*x_incr_t::del_j(),
            x_init_t::del_k() + index*x_incr_t::del_k()>;
            
        using added_sten_t = stencil_t<accum_t::center(), elem_t<xi_t, infos_t>>;
        using uni_t        = stencil_union<accum_t, added_sten_t>;
        using type = pattern_impl_t<x_init_t, x_incr_t, remaining-1, num_points, infos_t, uni_t>::type;
    };
    
    template <typename x_init_t, typename x_incr_t, const int num_points, typename infos_t, typename accum_t>
    requires (num_points > 0)
    struct pattern_impl_t<x_init_t, x_incr_t, 0, num_points, infos_t, accum_t>
    {
        using type = accum_t;
    };
    
    
    
    template <const grid::array_centering ctr, typename x_init_t, typename x_incr_t, const int num_points, typename infos_t>
    using pattern_t = typename pattern_impl_t<x_init_t, x_incr_t, num_points, num_points, infos_t, stencil_t<ctr>>::type;
}