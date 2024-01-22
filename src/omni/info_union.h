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

    // Takes the union of info lists
    template <typename... infos_t>
    using info_union = typename detail::multi_union_t<infos_t...>::type;
    
    namespace detail
    {
        template <typename kernel_t, typename... kernels_t>
        struct combine_info_impl
        {
            using type = info_union<typename kernel_t::info_type, typename combine_info_impl<kernels_t...>::type>;
        };

        template <typename kernel_t>
        struct combine_info_impl<kernel_t>
        {
            using type = typename kernel_t::info_type;
        };
    }
    
    namespace detail
    {
        template <typename cur_list_t, typename... infos_t>
        struct shmem_only;
        
        template <typename cur_list_t, typename info_t, typename... infos_t>
        struct shmem_only<cur_list_t, info_t, infos_t...>
        {
            using loc_type = typename std::conditional<info_t::is_shmem_buffered, typename cur_list_t::expand_list<info_t>, cur_list_t>::type;
            using type     = typename shmem_only<loc_type, infos_t...>::type;
        };
        
        template <typename cur_list_t, typename info_t>
        struct shmem_only<cur_list_t, info_t>
        {
            using loc_type = typename std::conditional<info_t::is_shmem_buffered, typename cur_list_t::expand_list<info_t>, cur_list_t>::type;
            using type     = loc_type;
        };
        
        template <typename nothing_t>
        struct shmem_info_impl;
        
        template <>
        struct shmem_info_impl<info_list_t<>>
        {
            using type = info_list_t<>;
        };
        
        template <typename... infos_t>
        struct shmem_info_impl<info_list_t<infos_t...>>
        {
            using type = shmem_only<info_list_t<>, infos_t...>::type;
        };
    }

    // Combined the "info type" of a series of types
    template <typename... kernels_t> using combine_infos = typename detail::combine_info_impl<kernels_t...>::type;
    
    // Given an info list, extracts only the shmem-bufferable infos
    template <typename list_t>
    using shmem_info_list = typename detail::shmem_info_impl<list_t>::type;
    
    
}