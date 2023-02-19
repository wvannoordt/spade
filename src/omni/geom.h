#pragma once

#include "core/grid_index_types.h"

namespace spade::omni
{
    template <const int di, const int dj, const int dk> struct offset_t
    {
        constexpr static bool is_1d = ((dj==0) && (di==0));
        
        //if we are at a node of type 'center', what kind of node is at the given offset?
        template <const grid::array_centering center>
        constexpr static grid::array_centering relative_node_centering = grid::cell_centered;
    };
}