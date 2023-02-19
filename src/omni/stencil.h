#pragma once

#include "core/grid_index_types.h"
#include "core/utils.h"

#include "omni/info.h"
#include "omni/geom.h"

namespace spade::omni
{
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
        
        template <const int idx> using stencil_element = typename utils::get_pack_type<idx, idx_list_t...>::type;
        constexpr static int num_elements() {return sizeof...(idx_list_t);}
        
        constexpr static grid::array_centering center() {return center_val;}
        template <const grid::array_centering ctr> static constexpr int count_type = 0;
        constexpr static int num_face() {return count_type<grid::face_centered>;}
        constexpr static int num_cell() {return count_type<grid::cell_centered>;}
        constexpr static int num_node() {return count_type<grid::node_centered>;}
    };
}