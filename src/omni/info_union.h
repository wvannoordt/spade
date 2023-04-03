#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace detail
    {
        template <typename info1_t, typename info2_t, const int elems_left>
        struct info_union_t
        {
            constexpr static int idx  = info2_t::num_infos() - elems_left;
            using cur_t               = typename info2_t::template info_elem<idx>;
            using true_t              = typename info_union_t<info1_t, info2_t, elems_left-1>::type;
            using false_t             = typename info_union_t<typename info1_t::template expand_list<cur_t>, info2_t, elems_left-1>::type;
            constexpr static bool add = info1_t::template contains<cur_t>;
            using type                = std::conditional<add, true_t, false_t>::type;
        };

        //Note: the trick with the '0' in the end of the partial specialization is nice for
        //sidestepping those silly CUDA compiler errors. Also prevents specialization for empty lists.
        template <typename info1_t, typename info2_t> struct info_union_t<info1_t, info2_t, 0>
        {
            using type = info1_t;
        };

        //union of two info lists
        template <typename info1_t, typename info2_t>
        using info_binary_union = typename info_union_t<info1_t, info2_t, info2_t::num_infos()>::type;

        template <typename info_t, typename... infos_t> struct multi_union_t
        {
            using type = info_binary_union<info_t, typename multi_union_t<infos_t...>::type>;
        };

        template <typename info1_t, typename info2_t> struct multi_union_t<info1_t, info2_t>
        {
            using type = info_binary_union<info1_t, info2_t>;
        };
    }

    //Only useful thing in this header
    template <typename... infos_t>
    using info_union = typename detail::multi_union_t<infos_t...>::type;
}