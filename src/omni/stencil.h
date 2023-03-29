#pragma once

#include "core/grid_index_types.h"
#include "core/utils.h"

#include "omni/info.h"
#include "omni/geom.h"

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
            const int idx,
            const int ar_size,
            const grid::array_centering ctr,
            const grid::array_centering query,
            typename idx_elem_t,
            typename... idx_list_t
        >
        struct node_count_helper_t
        {
            constexpr static int value = icount_t<ctr, query, idx_elem_t>::value + node_count_helper_t<idx+1, ar_size, ctr, query, idx_list_t...>::value;
        };
        
        template
        <
            const int ii,
            const grid::array_centering ctr,
            const grid::array_centering query,
            typename idx_elem_t
        >
        struct node_count_helper_t<ii-1, ii, ctr, query, idx_elem_t>
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
        static_assert(sizeof...(idx_list_t)>0, "Cannot specify a stencil with no nodes!");
        constexpr static bool is_directed()
        {
            // return idx_list_t::offset::is_1d && ...;
            return false;
        }
        
        template <const int idx> using stencil_element = typename utils::get_pack_type<idx, idx_list_t...>::type;
        constexpr static int num_elements() {return sizeof...(idx_list_t);}
        
        constexpr static grid::array_centering center() { return center_val; }
        
        template <const grid::array_centering ctr>
        static constexpr int count_type = detail::node_count_helper_t<0, sizeof...(idx_list_t), center_val, ctr, idx_list_t...>::value;
        
        constexpr static int num_face() {return count_type<grid::face_centered>;}
        constexpr static int num_cell() {return count_type<grid::cell_centered>;}
        constexpr static int num_node() {return count_type<grid::node_centered>;}
        constexpr static int num_edge() {return count_type<grid::edge_centered>;}
    };
}