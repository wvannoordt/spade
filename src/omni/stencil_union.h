#pragma once

#include "omni/stencil.h"
#include "omni/convert.h"

namespace spade::omni
{
    namespace detail
    {
        template <typename stencil1_t, typename stencil2_t, const int remaining> struct stencil_union_t
        {
            constexpr static bool st_eq = (stencil1_t::center() == stencil2_t::center()) && (stencil2_t::center() != grid::agno_centered);
            constexpr static bool ag_1  = (stencil1_t::center() == grid::agno_centered)  && (stencil2_t::center() != grid::agno_centered);
            constexpr static bool ag_2  = (stencil2_t::center() == grid::agno_centered)  && (stencil1_t::center() != grid::agno_centered);
            static_assert(st_eq || ag_1 || ag_2, "Attempted to union two incompatible stencils");

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
    
        template <const grid::array_centering c0, const grid::array_centering c1>
        struct most_specific_centering
        {
            //This is wrong, but the static assert in the stencil_union_t will give a better error message
            constexpr static grid::array_centering value = c0;
        };
        
        template <const grid::array_centering c0>
        requires (c0 != grid::agno_centered)
        struct most_specific_centering<c0, grid::agno_centered>
        {
            constexpr static grid::array_centering value = c0;
        };

        template <const grid::array_centering c0>
        requires (c0 != grid::agno_centered)
        struct most_specific_centering<grid::agno_centered, c0>
        {
            constexpr static grid::array_centering value = c0;
        };

        template <typename stencil1_t, typename stencil2_t>   using stencil_binary_union
        = typename detail::stencil_union_t<
            typename stencil1_t::recenter<most_specific_centering<stencil1_t::center(), stencil2_t::center()>::value>,
            stencil2_t,
            stencil2_t::num_elements()
            >::type;

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

    namespace detail
    {
        template <typename kernel_t, typename... kernels_t>
        struct combine_omni_impl
        {
            using type = stencil_union<typename kernel_t::omni_type, typename combine_omni_impl<kernels_t...>::type>;
        };

        template <typename kernel_t>
        struct combine_omni_impl<kernel_t>
        {
            using type = typename kernel_t::omni_type;
        };
    }

    template <typename... kernels_t> using combine_omni_stencils = typename detail::combine_omni_impl<kernels_t...>::type;
}