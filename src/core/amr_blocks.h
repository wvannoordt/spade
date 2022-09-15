#pragma once

#include <vector>

#include "core/ctrs.h"
#include "core/amr_node.h"
#include "core/block_config.h"

namespace spade::block_config
{
    template <typename coord_val_t, const std::size_t grid_dim, const amr::amr_block_count_mode count_mode = amr::amr_count_terminal>
    struct amr_blocks_t
    {
        using node_type      = amr::amr_node_t<grid_dim>;
        using coord_val_type = coord_val_t;
        
        ctrs::array<int, grid_dim>                      num_blocks;
        bound_box_t<coord_val_t, grid_dim>              bounds;
        std::vector<node_type>                          root_nodes;
        std::vector<bound_box_t<coord_val_t, grid_dim>> block_boxes;
        std::vector<node_type*>                         all_nodes;
        
        amr_blocks_t(const ctrs::array<int, grid_dim>& num_blocks_in,
            const bound_box_t<coord_val_t, grid_dim>& bounds_in)
        {            
            num_blocks = num_blocks_in;
            bounds     = bounds_in;
            std::size_t nodes_size = 1;
            for (auto i: range(0, grid_dim)) nodes_size *= num_blocks[i];
            root_nodes.resize(nodes_size);
            for (auto i: range(0, root_nodes.size()))
            {
                auto ijk = ctrs::expand_index(i, num_blocks);
                auto num_periphs   = static_math::pow<3,grid_dim>::value;
                auto num_neighbors = num_periphs - 1;
                root_nodes[i].neighbors.reserve(num_neighbors);
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
                        neighbor_relationship.endpoint = &root_nodes[j];
                        neighbor_relationship.edge     = dijk;
                        root_nodes[i].neighbors.push_back(neighbor_relationship);
                    }
                }
                root_nodes[i].parent = nullptr;
                for (auto d: range(0, grid_dim))
                {
                    root_nodes[i].amr_position.min(d) = amr::amr_coord_t(ijk[d],   num_blocks[d], 0);
                    root_nodes[i].amr_position.max(d) = amr::amr_coord_t(ijk[d]+1, num_blocks[d], 0);
                }
            }
            this->enumerate();
        }
        
        //Note: these return the block bounding box in computational coordinates
        bound_box_t<coord_val_t, grid_dim>& get_block_box (const std::size_t& lb_glob)
        {
            return block_boxes[lb_glob];
        }
        
        //Note: these return the block bounding box in computational coordinates
        const bound_box_t<coord_val_t, grid_dim>& get_block_box (const std::size_t& lb_glob) const
        {
            return block_boxes[lb_glob];
        }
        
        std::size_t get_num_blocks() const { return all_nodes.size(); }
        template <typename filter_t> std::size_t get_num_blocks(const filter_t& filter) const
        {
            std::size_t output = 0;
            for (const auto p: all_nodes)
            {
                if (filter(*p)) ++output;
            }
            return output;
        }
        
        static constexpr std::size_t dim() { return grid_dim; }
        
        // re-collects the global list of amr nodes and their bounding boxes,
        // must be called after every single refinement operation
        void enumerate()
        {
            std::size_t block_count = 0;
            for (auto& n: root_nodes)
            {
                block_count += n.template count_nodes<count_mode>();
            }
            all_nodes.reserve(block_count);
            for (auto& n: root_nodes)
            {
                n.template collect_nodes<count_mode>(all_nodes);
            }
            block_boxes.resize(block_count);
            ctrs::array<coord_val_t, grid_dim> dx;
            for (auto i: range(0, grid_dim)) dx[i] = bounds.size(i)/num_blocks[i];
            for (auto i: range(0, block_boxes.size()))
            {
                const auto& amr_bbox  = all_nodes[i]->amr_position;
                auto& comp_bbox = block_boxes[i];
                for (auto d: range(0, grid_dim))
                {
                    comp_bbox.min(d) = amr_bbox.min(d).convert_to_coordinate(bounds.min(d), dx[d]);
                    comp_bbox.max(d) = amr_bbox.max(d).convert_to_coordinate(bounds.min(d), dx[d]);
                }
            }
        }
    };
}