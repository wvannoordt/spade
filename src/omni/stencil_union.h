#pragma once

#include "omni/stencil.h"

namespace spade::omni
{
    namespace detail
    {
        template <typename stencil1_t, typename stencil2_t, const int remaining> struct stencil_union_t
        {
            constexpr static int index = stencil2_t::num_elements() - remaining;
            using element_t            = stencil2_t::template stencil_element<index>;
            using query_offset_t       = typename element_t::position_type;
            using query_info_t         = typename element_t::info_type;

            constexpr static bool s1_has_element = stencil1_t::template contains_at<query_offset_t>;
            using true_t                         = stencil_insert_at<stencil1_t, query_offset_t, query_info_t>;
            using false_t                        = stencil1_t::template extend<elem_t<query_offset_t, query_info_t>>;
            using merged_t                       = typename std::conditional<s1_has_element, true_t, false_t>::type;
            using type                           = typename stencil_union_t<merged_t, stencil2_t, remaining-1>::type;
        };

        template <typename stencil1_t, typename stencil2_t> struct stencil_union_t<stencil1_t, stencil2_t, 0>
        {
            using type = stencil1_t;
        };

        template <typename stencil1_t, typename stencil2_t>   using stencil_binary_union = typename detail::stencil_union_t<stencil1_t, stencil2_t, stencil2_t::num_elements()>::type;

        template <typename stencil_t, typename... stencils_t> struct multi_sten_union_t
        {
            using type = stencil_binary_union<stencil_t, typename multi_sten_union_t<stencils_t...>::type>;
        };
        template <typename stencil_t> struct multi_sten_union_t<stencil_t>
        {
            using type = stencil_t;
        };
    }

    template <typename... stencils_t> using stencil_union = typename detail::multi_sten_union_t<stencils_t...>::type;
}