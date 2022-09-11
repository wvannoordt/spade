#pragma once

#include <vector>

#include "core/bounding_box.h"
#include "core/amr_coord.h"

namespace spade::amr
{
    template <const std::size_t grid_dim> struct amr_node_t
    {
        bound_box_t<amr_coord_t, grid_dim> amr_position;
        std::vector<amr_node_t> subnodes;
        amr_node_t(){}
    };
}