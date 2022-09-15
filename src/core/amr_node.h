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
        bool terminal() const { return subnodes.size()==0; }
        
        template <const amr_block_count_mode count_mode>
        std::size_t count_nodes() const
        {
            if constexpr (count_mode==amr_count_terminal)
            {
                if (this->terminal()) return 1;
            }
            std::size_t output = 0;
            if constexpr (count_mode==amr_count_all) ++output;
            for (auto& n: subnodes) output += n.template count_nodes<count_mode>();
            return output;
        }
        
        template <const amr_block_count_mode count_mode>
        void collect_nodes(std::vector<amr_node_t*>& node_list)
        {
            node_list.push_back(this);
            for (auto& n: subnodes) n.template collect_nodes<count_mode>(node_list);
        }
        
        amr_node_t(){}
    };
}