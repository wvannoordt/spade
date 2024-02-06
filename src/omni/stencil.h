#pragma once

#include "grid/grid_index_types.h"
#include "core/utils.h"

#include "omni/info.h"
#include "omni/geom.h"
#include "omni/info_union.h"

namespace spade::omni
{
    namespace detail
    {
        template <const grid::array_centering ctr, const grid::array_centering query, typename idx_elem_t>
        struct icount_t
        {
            constexpr static int value = (query==idx_elem_t::position_type::template relative_node_centering<ctr>)?1:0;
        };
        
        template
        <
            const grid::array_centering ctr,
            const grid::array_centering query,
            typename idx_elem_t,
            typename... idx_list_t
        >
        struct node_count_helper_t
        {
            constexpr static int value = icount_t<ctr, query, idx_elem_t>::value + node_count_helper_t<ctr, query, idx_list_t...>::value;
        };
        
        template
        <
            const grid::array_centering ctr,
            const grid::array_centering query,
            typename idx_elem_t
        >
        struct node_count_helper_t<ctr, query, idx_elem_t>
        {
            constexpr static int value = icount_t<ctr, query, idx_elem_t>::value;
        };
    }
    
    template <typename position_t, typename info_t> struct elem_t
    {
        using position_type = position_t;
        using info_type     = info_t;
    };

    namespace detail
    {
        template <typename offset_query_t, typename... info_lists_t>
        struct search_info_impl;
        
        template <typename offset_query_t, typename info_list_t, typename... info_lists_t>
        struct search_info_impl<offset_query_t, info_list_t, info_lists_t...>
        {
            using type = typename std::conditional<
                std::same_as<typename info_list_t::position_type, offset_query_t>,
                typename info_list_t::info_type,
                typename search_info_impl<offset_query_t, info_lists_t...>::type>::type;
        };
        
        template <typename offset_query_t, typename info_list_t>
        struct search_info_impl<offset_query_t, info_list_t>
        {
            using type = typename info_list_t::info_type;
        };
        
        
        template <const int multi, const int val, const int new_val>
        struct extent_reduce_impl
        {
            constexpr static int value = (multi*new_val > multi*val) ? new_val : val;
        };
        
        template <const int multi, const int idx, const int current, typename... elems_t>
        struct extent_impl;
        
        template <const int multi, const int idx, const int current, typename offset_t, typename list_t, typename... elems_t>
        struct extent_impl<multi, idx, current, elem_t<offset_t, list_t>, elems_t...>
        {
            constexpr static int nval        = offset_t::template elem<idx>;
            constexpr static int new_current = extent_reduce_impl<multi, current, nval>::value;
            constexpr static int value       = extent_impl<multi, idx, new_current, elems_t...>::value;
        };
        
        template <const int multi, const int idx, const int current, typename offset_t, typename list_t>
        struct extent_impl<multi, idx, current, elem_t<offset_t, list_t>>
        {
            constexpr static int value = extent_reduce_impl<multi, current, offset_t::template elem<idx>>::value;
        };
    }
    
    template <const grid::array_centering center_val, typename... idx_list_t>
    struct stencil_t
    {
        template <typename offset_query_t> static constexpr bool contains_at = (std::same_as<typename idx_list_t::position_type, offset_query_t> || ...);
        
        template <typename offset_query_t> using info_at
            = typename std::conditional<contains_at<offset_query_t>, typename detail::search_info_impl<offset_query_t, idx_list_t...>::type, info_list_t<>>::type;
        
        template <const int idx> using stencil_element = typename utils::get_pack_type<idx, idx_list_t...>::type;
        constexpr static int num_elements() {return sizeof...(idx_list_t);}
        
        constexpr static grid::array_centering center() { return center_val; }
        
        template <const grid::array_centering ctr>
        static constexpr int count_type = detail::node_count_helper_t<center_val, ctr, idx_list_t...>::value;
        
        constexpr static int num_face() {return count_type<grid::face_centered>;}
        constexpr static int num_cell() {return count_type<grid::cell_centered>;}
        constexpr static int num_node() {return count_type<grid::node_centered>;}
        constexpr static int num_edge() {return count_type<grid::edge_centered>;}

        //adds a stencil point to the list (even if it already exists)
        template <typename new_elem_t> using extend = stencil_t<center_val, idx_list_t..., new_elem_t>;
        template <const grid::array_centering new_center> using recenter = stencil_t<new_center, idx_list_t...>;
        
        template <const int multi, const int idx> constexpr static int extent_impl = detail::extent_impl<multi, idx, -100*multi, idx_list_t...>::value;
        
        template <const int idx> constexpr static int max_extent = extent_impl< 1, idx>;
        template <const int idx> constexpr static int min_extent = extent_impl<-1, idx>;
    };

    namespace detail
    {
        template <typename sten_t, typename off_t, typename list_t, const int remaining, typename accum_t>
        struct sten_insert_list_impl
        {
            constexpr static int index     = sten_t::num_elements() - remaining;
            using query_elem_t             = sten_t::template stencil_element<index>;
            using cur_offset_t             = query_elem_t::position_type;
            using cur_list_t               = query_elem_t::info_type;
            constexpr static bool is_query = std::same_as<off_t, cur_offset_t>;
            using new_list_t               = typename std::conditional<is_query, info_union<cur_list_t, list_t>, cur_list_t>::type;
            using new_elem_t               = elem_t<cur_offset_t, new_list_t>;
            using new_accum_t              = accum_t::template extend<new_elem_t>;
            using type                     = typename sten_insert_list_impl<sten_t, off_t, list_t, remaining - 1, new_accum_t>::type;
        };

        template <typename sten_t, typename off_t, typename list_t, typename accum_t>
        struct sten_insert_list_impl<sten_t, off_t, list_t, 0, accum_t>
        {
            using type = accum_t;
        };

        //There is almost certainly a better way of doing this
        template <typename istencil_t, const grid::array_centering ctr, const int idx, const int found, const int remaining>
        struct stencil_seek_offset_impl
        {
            static_assert(idx < istencil_t::template count_type<ctr>, "seek index must not exceed the number of elements on stencil");
            constexpr static int query_idx = istencil_t::num_elements() - remaining;
            //if the centering of the query element matches the centering being sought, but not right index
            using match_type   = typename stencil_seek_offset_impl<istencil_t, ctr, idx, found+1, remaining-1>::type;

            //if the centering of the query element doesn't match the centering being sought
            using nomatch_type = typename stencil_seek_offset_impl<istencil_t, ctr, idx, found,   remaining-1>::type;

            using element_type                = istencil_t::template stencil_element<query_idx>;
            using offset_type                 = typename element_type::position_type;

            constexpr static grid::array_centering query_centering 
                = offset_type::template relative_node_centering<istencil_t::center()>;

            constexpr static bool center_match = (ctr == query_centering);
            constexpr static bool total_match  = center_match && (found == idx);

            using next_type = typename std::conditional<center_match, match_type, nomatch_type>::type;
            using end_type  = offset_type;
            using type = typename std::conditional<total_match, end_type, next_type>::type;
        };

        template <typename istencil_t, const grid::array_centering ctr, const int idx, const int found>
        struct stencil_seek_offset_impl<istencil_t, ctr, idx, found, 0>
        {
            //invalid, never realized thanks to the static_assert in the general case
            using type = int;
        };
    }

    template <typename istencil_t, typename cur_offset_t, typename new_info_list_t> using stencil_insert_at
        = typename detail::sten_insert_list_impl<istencil_t, cur_offset_t, new_info_list_t, istencil_t::num_elements(), stencil_t<istencil_t::center()>>::type;

    template <typename istencil_t, const grid::array_centering ctr, const int idx> using offset_at
        = typename detail::stencil_seek_offset_impl<istencil_t, ctr, idx, 0, istencil_t::num_elements()>::type;
    
    namespace detail
    {
        // This is one of the worst implementations I have ever constructed. If the compiler melts a machine,
        // It could very well be this implementation
        template <typename istencil_t, const grid::array_centering ctr, typename offset_query_t, const int tsize, const int remaining>
        struct stencil_index_of_offset_impl
        {
            constexpr static int idx    = istencil_t::template count_type<ctr> - remaining;
            using elem_offset_type      = offset_at<istencil_t, ctr, idx>;
            constexpr static bool found = std::same_as<offset_query_t, elem_offset_type>;
            constexpr static int value  = found ? idx : stencil_index_of_offset_impl<istencil_t, ctr, offset_query_t, tsize, remaining-1>::value;
        };

        template <typename istencil_t, const grid::array_centering ctr, typename offset_query_t, const int tsize>
        struct stencil_index_of_offset_impl<istencil_t, ctr, offset_query_t, tsize, 0>
        {
            constexpr static int value  = -1;
        };
    }

    template <typename istencil_t, typename offset_query_t> constexpr static int index_of
        = detail::stencil_index_of_offset_impl<
            istencil_t,
            offset_query_t::template relative_node_centering<istencil_t::center()>,
            offset_query_t,
            istencil_t::template count_type<offset_query_t::template relative_node_centering<istencil_t::center()>>,
            istencil_t::template count_type<offset_query_t::template relative_node_centering<istencil_t::center()>>
            >::value;
}