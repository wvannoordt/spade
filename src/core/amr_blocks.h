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
            for (auto i: range(0, nodes.size()))
            {
                auto ijk = ctrs::expand_index(i, num_blocks);
                auto num_periphs   = static_math::pow<3,grid_dim>::value;
                auto num_neighbors = num_periphs - 1;
                nodes[i].neighbors.reserve(num_neighbors);
                for (auto jt: range(0,num_periphs))
                {
                    ctrs::array<int,grid_dim> extent = 3;
                    auto dijk = ctrs::expand_index(jt, extent);
                    for (auto& e: dijk) e -= 1;
                    bool all_zero = true;
                    for (const auto& e: dijk) all_zero = (all_zero && (e==0));
                    if (!all_zero)
                    {
                        auto ijk_neigh = dijk + ijk;
                        for (auto k: range(0, grid_dim))
                        {
                            if (ijk_neigh[k]<0)              ijk_neigh[k] += num_blocks[k];
                            if (ijk_neigh[k]>=num_blocks[k]) ijk_neigh[k] -= num_blocks[k];
                        }
                        auto j = collapse_index(ijk_neigh, num_blocks);
                        
                        //node i has neighbor j
                        amr::amr_neighbor_t<grid_dim> neighbor_relationship;
                        neighbor_relationship.endpoint = &nodes[j];
                        neighbor_relationship.edge     = dijk;
                        nodes[i].neighbors.push_back(neighbor_relationship);
                    }
                }
                nodes[i].parent = nullptr;
                for (auto d: range(0, grid_dim))
                {
                    nodes[i].amr_position.min(d) = amr::amr_coord_t(ijk[d],                   0);
                    nodes[i].amr_position.max(d) = amr::amr_coord_t((ijk[d]+1)%num_blocks[d], 0);
                }
            }
        }
        
        void enumerate()
        {
            
        }
    };
}