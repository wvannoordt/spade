#pragma once

#include <vector>

#include "core/ctrs.h"
#include "core/bounding_box.h"
#include "core/amr_coord.h"

namespace spade::amr
{
    enum amr_block_count_mode
    {
        amr_count_terminal,
        amr_count_all
    };
    
    template <const std::size_t grid_dim> struct amr_node_t;
    template <const std::size_t grid_dim> struct amr_neighbor_t
    {
        //                            edge
        //         amr_node          ------>     endpoint
        amr_node_t<grid_dim>* endpoint;
        ctrs::array<int, grid_dim> edge;
        
    };
    template <const std::size_t grid_dim> struct amr_node_t
    {
        bound_box_t<amr_coord_t, grid_dim> amr_position;
        std::vector<amr_node_t> subnodes;
        std::vector<amr_neighbor_t<grid_dim>> neighbors;
        amr_node_t<grid_dim>* parent;
        amr_node_t(){}
    };
}