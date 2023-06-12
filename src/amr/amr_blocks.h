#pragma once

#include <vector>

#include "core/ctrs.h"
#include "amr/amr_node.h"
#include "core/block_config.h"

namespace spade::amr
{
    template <typename coord_val_t, typename array_designator_t>
    struct amr_blocks_t
    {
        using node_type      = amr::amr_node_t<array_designator_t::size()>;
        using coord_val_type = coord_val_t;
        
        ctrs::array<int, 3>                      num_blocks;
        bound_box_t<coord_val_t, 3>              bounds;
        std::vector<node_type>                   root_nodes;
        std::vector<bound_box_t<coord_val_t, 3>> block_boxes;
        std::vector<node_type*>                  all_nodes;
        
        constexpr static std::size_t dim() { return array_designator_t::size(); }
        
        amr_blocks_t(const array_designator_t& num_blocks_in,
            const bound_box_t<coord_val_t, array_designator_t::size()>& bounds_in)
        {
            ctrs::copy_array(num_blocks_in, num_blocks, 1);
            bounds.min(2) = 0.0;
            bounds.max(2) = 1.0;
            for (auto d: range(0, dim()))
            {
                bounds.min(d) = bounds_in.min(d);
                bounds.max(d) = bounds_in.max(d);
            }
            std::size_t nodes_size = 1;
            for (auto i: range(0, dim())) nodes_size *= num_blocks[i];
            root_nodes.resize(nodes_size);
            for (auto i: range(0, root_nodes.size()))
            {
                auto ijk = ctrs::expand_index(i, num_blocks);
                auto num_periphs   = static_math::pow<3,dim()>::value;
                auto num_neighbors = num_periphs - 1;
                root_nodes[i].neighbors.reserve(num_neighbors);
                for (auto jt: range(0,num_periphs))
                {
                    ctrs::array<int,3> extent = 3;
                    if constexpr (dim() == 2) extent[2] = 1;
                    auto dijk = ctrs::expand_index(jt, extent);
                    for (auto& e: dijk) e -= 1;
                    bool all_zero = true;
                    for (const auto& e: dijk) all_zero = (all_zero && (e==0));
                    if (!all_zero)
                    {
                        auto ijk_neigh = dijk + ijk;
                        for (auto k: range(0, dim()))
                        {
                            if (ijk_neigh[k]<0)              ijk_neigh[k] += num_blocks[k];
                            if (ijk_neigh[k]>=num_blocks[k]) ijk_neigh[k] -= num_blocks[k];
                        }
                        auto j = collapse_index(ijk_neigh, num_blocks);
                        
                        //node i has neighbor j
                        amr::amr_neighbor_t<dim()> neighbor_relationship;
                        neighbor_relationship.endpoint = &root_nodes[j];
                        neighbor_relationship.edge     = dijk;
                        root_nodes[i].neighbors.push_back(neighbor_relationship);
                    }
                }
                root_nodes[i].parent = nullptr;
                for (auto d: range(0, dim()))
                {
                    root_nodes[i].amr_position.min(d) = amr::amr_coord_t(ijk[d],   num_blocks[d], 0);
                    root_nodes[i].amr_position.max(d) = amr::amr_coord_t(ijk[d]+1, num_blocks[d], 0);
                }
                root_nodes[i].level = 0;
            }
            this->enumerate();
        }
        
        //Note: these return the block bounding box in computational coordinates
        const auto& get_block_box (const std::size_t& lb_glob) const
        {
            return block_boxes[lb_glob];
        }
        
        std::size_t total_num_blocks() const { return all_nodes.size(); }
        template <typename filter_t> std::size_t get_num_blocks(const filter_t& filter) const
        {
            std::size_t output = 0;
            for (const auto p: all_nodes)
            {
                if (filter(*p)) ++output;
            }
            return output;
        }
        
        // re-collects the global list of amr nodes and their bounding boxes,
        // must be called after every single refinement operation
        void enumerate()
        {
            std::size_t block_count = 0;
            for (auto& n: root_nodes)
            {
                block_count += n.count_nodes();
            }
            all_nodes.clear();
            all_nodes.reserve(block_count);
            for (auto& n: root_nodes)
            {
                n.collect_nodes(all_nodes);
            }
            block_boxes.clear();
            block_boxes.resize(block_count);
            ctrs::array<coord_val_t, dim()> dx;
            for (auto i: range(0, dim())) dx[i] = bounds.size(i)/num_blocks[i];
            for (auto i: range(0, block_boxes.size()))
            {
                const auto& amr_bbox  = all_nodes[i]->amr_position;
                auto& comp_bbox = block_boxes[i];
                for (auto d: range(0, dim()))
                {
                    comp_bbox.min(d) = amr_bbox.min(d).convert_to_coordinate(bounds.min(d), dx[d]);
                    comp_bbox.max(d) = amr_bbox.max(d).convert_to_coordinate(bounds.min(d), dx[d]);
                }
            }
        }
        
        const auto&       get_bounding_box(const std::size_t lb)              const { return block_boxes[lb]; }
        const coord_val_t get_size(const std::size_t i, const std::size_t lb) const { return block_boxes[lb].size(i); }
    };
}