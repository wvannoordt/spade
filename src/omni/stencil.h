#pragma once

#include "core/grid_index_types.h"
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

    template <const grid::array_centering center_val, typename... idx_list_t>
    struct stencil_t
    {
        constexpr static bool is_directed()
        {
            // return idx_list_t::offset::is_1d && ...;
            return false;
        }
        template <typename offset_query_t> static constexpr bool contains_at = (std::same_as<typename idx_list_t::position_type, offset_query_t> || ...);
        
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
    }

    template <typename istencil_t, typename cur_offset_t, typename new_info_list_t> using stencil_insert_at
        = typename detail::sten_insert_list_impl<istencil_t, cur_offset_t, new_info_list_t, istencil_t::num_elements(), stencil_t<istencil_t::center()>>::type;
}