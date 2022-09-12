#pragma once

#include <vector>

#include "core/ctrs.h"
#include "core/amr_node.h"
#include "core/block_config.h"

namespace spade::block_config
{
    template <typename coord_val_t, const std::size_t grid_dim>
    struct amr_blocks_t
    {
        using node_type = amr::amr_node_t<grid_dim>;
        ctrs::array<int, grid_dim> num_blocks;
        bound_box_t<coord_val_t, grid_dim>  bounds;
        std::vector<node_type> nodes;
        
        amr_blocks_t(const ctrs::array<int, grid_dim>& num_blocks_in,
            const bound_box_t<coord_val_t, grid_dim>& bounds_in)
        {
            num_blocks = num_blocks_in;
            bounds     = bounds_in;
            std::size_t nodes_size = 1;
            for (auto i: range(0, grid_dim)) nodes_size *= num_blocks[i];
            nodes.resize(nodes_size);
        }
        
        void enumerate()
        {
            
        }
    };
}