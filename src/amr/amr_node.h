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
        using amr_refine_t = ctrs::array<bool,grid_dim>;
        bound_box_t<amr_coord_t, grid_dim> amr_position;
        std::vector<amr_node_t> subnodes;
        std::vector<amr_neighbor_t<grid_dim>> neighbors;
        ctrs::array<int, grid_dim> level;
        amr_node_t<grid_dim>* parent;
        
        static constexpr std::size_t dim() { return grid_dim; }
        
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
            
            for (int i = 0; i < subnodes.size(); ++i)
            {
                auto& child = subnodes[i];
                
                //First add neighbor relationships for siblings of current child
                for (int j = 0; j < subnodes.size(); ++j)
                {
                    if (i != j)
                    {
                        auto& sibling = subnodes[j];
                        child.create_neighbor(sibling);
                    }
                }
                
                //Now add parent neighbors
                for (const auto& neigh: neighbors)
                {
                    child.create_neighbor(*neigh.endpoint);
                }
            }
        }
        
        void create_neighbor(amr_node_t& candidate)
        {
            if (!candidate.terminal()) return;
            ctrs::array<ctrs::array<bool, 3>, 3> table;
            table.fill(true);
            for (int d = 0; d < dim(); ++d)
            {
                auto& in_region = table[d];
                const auto& n_min = this->amr_position.min(d);
                const auto& n_max = this->amr_position.max(d);
                
                const auto& m_min = candidate.amr_position.min(d);
                const auto& m_max = candidate.amr_position.max(d);
                
                in_region[0] = (m_max >= n_min) && (m_min <  n_min); // this one
                in_region[1] = (m_max >  n_min) && (m_min <  n_max);
                in_region[2] = (m_max >  n_max) && (m_min <= n_max); // and this one
            }
            
            int kmax = (dim()==3)?3:1;
            for (int k = 0; k < kmax; ++k)
            {
                for (int j = 0; j < 3; ++j)
                {
                    for (int i = 0; i < 3; ++i)
                    {
                        ctrs::array<int, 3> edge(i-1, j-1, k-1);
                        if constexpr (dim() == 2) edge[2] = 0;
                        
                        bool valid = table[0][i] && table[1][j] && table[2][k];
                        if (valid)
                        {
                            amr_neighbor_t<dim()> neighbor;
                            neighbor.endpoint = &candidate;
                            ctrs::copy_array(edge, neighbor.edge);
                            neighbors.push_back(neighbor);
                        }
                    }
                }
            }
        }
        
        void collect_nodes(std::vector<amr_node_t*>& node_list)
        {
            node_list.push_back(this);
            for (auto& n: subnodes) n.collect_nodes(node_list);
        }
        
        auto ubox() const
        {
            bound_box_t<double, grid_dim> out;
            for (int d = 0; d < grid_dim; ++d)
            {
                out.min(d) = amr_position.min(d).to_unit();
                out.max(d) = amr_position.max(d).to_unit();
            }
            return out;
        }
        
        bound_box_t<bool, 3> is_domain_boundary() const
        {
            bound_box_t<bool, 3> out;
            out.min(2) = true;
            out.max(2) = true;
            for (int i = 0; i < dim(); ++i)
            {
                out.min(i) = (amr_position.min(i).partition == 0) && (amr_position.min(i).bits == 0);
                out.max(i) = (amr_position.max(i).partition == amr_position.max(i).num_partitions);
            }
            return out;
        }
        
        void debug() const
        {
            print(neighbors.size(), this->ubox());
            for (const auto& i:neighbors)
            {
                print("----------------------");
                print(i.edge);
                print(i.endpoint->ubox());
                print("----------------------");
            }
        }
        
        amr_node_t(){}
    };
}