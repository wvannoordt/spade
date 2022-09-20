#pragma once

#include <vector>

#include "core/ctrs.h"
#include "core/bounding_box.h"
#include "amr/amr_coord.h"

namespace spade::amr
{    
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
        using amr_refine_t = ctrs::array<bool,3>;
        bound_box_t<amr_coord_t, grid_dim> amr_position;
        std::vector<amr_node_t> subnodes;
        std::vector<amr_neighbor_t<grid_dim>> neighbors;
        ctrs::array<int, grid_dim> level;
        amr_node_t<grid_dim>* parent;
        
        static constexpr std::size_t dim() {return grid_dim;}
        
        bool terminal() const { return subnodes.size()==0; }
        
        std::size_t count_nodes() const
        {
            std::size_t output = 1;
            for (auto& n: subnodes) output += n.count_nodes();
            return output;
        }
        
        void refine_node_recurse(const amr_refine_t& ref_type)
        {
            std::size_t num_sub_nodes = 1;
            for (auto i: range(0, dim())) num_sub_nodes *= ref_type[i]?2:1;
            //exception conditions
            if (num_sub_nodes == 1) return;
            if (!this->terminal()) return;
            subnodes.resize(num_sub_nodes);
            auto rg = range(0,ref_type[0]?2:1)*range(0,ref_type[1]?2:1)*range(0,(ref_type[dim()-1]&&(dim()==3))?2:1);
            int ct = 0;
            for (auto i: rg)
            {
                auto& child = subnodes[ct];
                child.parent = this;
                child.amr_position = this->amr_position;
                for (auto d: range(0, dim()))
                {
                    child.level[d] = this->level[d];
                    if (ref_type[d])
                    {
                        child.level[d]++;
                        amr_coord_t mid_coord = this->amr_position.min(d);
                        mid_coord.set_bit(child.level[d]-1, 1);
                        if (i[d]==1) child.amr_position.min(d) = mid_coord;
                        else         child.amr_position.max(d) = mid_coord;
                    }
                }
                ct++;
            }
        }
        
        void collect_nodes(std::vector<amr_node_t*>& node_list)
        {
            node_list.push_back(this);
            for (auto& n: subnodes) n.collect_nodes(node_list);
        }
        
        amr_node_t(){}
    };
}