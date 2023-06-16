#pragma once

#include <vector>
#include <tuple>

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
        amr_node_t<grid_dim>* endpoint = nullptr;
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
        bool flood_lock = false;
        
        static constexpr std::size_t dim() { return grid_dim; }
        
        bool terminal() const { return subnodes.size()==0; }
        
        std::size_t count_nodes() const
        {
            std::size_t output = 1;
            for (auto& n: subnodes) output += n.count_nodes();
            return output;
        }
        
        void refine_node_recurse(const amr_refine_t& ref_type, const ctrs::array<bool, dim()> is_periodic = false)
        {
            std::size_t num_sub_nodes = 1;
            for (auto i: range(0, dim())) num_sub_nodes *= ref_type[i]?2:1;
            
            //exception conditions
            if (num_sub_nodes == 1) return;
            if (!this->terminal())  return;
            if (this->flood_lock)   return;
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
                    auto cond =  [=](const auto& relationship) {return ctrs::dot_prod(neigh.edge, relationship.edge) >= 0;};
                    child.create_neighbor(*neigh.endpoint);
                    neigh.endpoint->create_neighbor(child);
                }
            }
            
            //Each of my neighbors track their own neighbors, so need to update those.
            //NOTE: the following two loops must be executed separately, but I am not quite sure why...
            //Remove any duplicates unduly created from periodicity
            for (auto& neigh: neighbors) neigh.endpoint->remove_duplicate_neighbors();
            
            //Now that I am no longer a terminal node, remove myself from my neighbors' list of neighbors
            for (auto& neigh: neighbors) neigh.endpoint->remove_neighbors([](const auto& ngh){ return !ngh.endpoint->terminal(); });
            
            //Generalize later: enforce no factor-4 refinement between neighboring nodes
            for (auto& child: subnodes) child.enforce_interface_compatibility(is_periodic);
        }
        
        void enforce_interface_compatibility(const ctrs::array<bool, dim()> is_periodic)
        {
            //todo: add in logic for derefinement
            this->flood_lock = true;
            const auto is_boundary = this->is_domain_boundary();
            const auto neighbor_copy = neighbors; //neighbors will change size during iteration otherwise
            for (auto& neigh: neighbor_copy)
            {
                //here we loop over neighbors to see if interface compatibility conditions are met
                amr_refine_t ref = false;
                bool refine_any  = false;
                for (int d = 0; d < dim(); ++d)
                {
                    if (this->level[d] > neigh.endpoint->level[d]+1) ref[d] = true;
                    refine_any = refine_any || ref[d];
                }
                
                //ignore if specified by the periodicity condition
                for (int d = 0; d < dim(); ++d)
                {
                    if (!is_periodic[d] && ((is_boundary.max(d) && neigh.edge[d] == 1) || (is_boundary.min(d) && neigh.edge[d] == -1))) refine_any = false;
                }
                if (refine_any)
                {
                    neigh.endpoint->refine_node_recurse(ref);
                }
            }
            this->flood_lock = false;
        }
        
        template <typename condition_t> void create_neighbor(amr_node_t& candidate, const condition_t& condition)
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
                
                const auto lboundary_wrap = [](const auto& x0, const auto& x1)
                {
                    auto y0 = x0;
                    auto y1 = x1;
                    if (y1.partition == y1.num_partitions)
                    {
                        y1.partition -= y1.num_partitions;
                        y0.partition -= y0.num_partitions;
                    }
                    return std::make_tuple(y0, y1);
                };
                
                const auto uboundary_wrap = [](const auto& x0, const auto& x1)
                {
                    auto y0 = x0;
                    auto y1 = x1;
                    if (y0.partition == 0 && y0.bits == 0)
                    {
                        y1.partition += y1.num_partitions;
                        y0.partition += y0.num_partitions;
                    }
                    return std::make_tuple(y0, y1);
                };
                
                //Wrap around boundaries to detect periodic neighbors
                const auto [m_min_p, m_max_p] = uboundary_wrap(m_min, m_max);
                const auto [m_min_m, m_max_m] = lboundary_wrap(m_min, m_max);
                
                in_region[0] = (m_max_m >= n_min) && (m_min_m <  n_min);
                in_region[1] = (m_max   >  n_min) && (m_min   <  n_max);
                in_region[2] = (m_max_p >  n_max) && (m_min_p <= n_max);
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
                            if (condition(neighbor)) neighbors.push_back(neighbor);
                        }
                    }
                }
            }
        }
        
        void create_neighbor(amr_node_t& candidate)
        {
            this->create_neighbor(candidate, [](const auto&i){return true;});
        }
        
        template <typename condition_t>
        void remove_neighbors(const condition_t& condition)
        {
            auto it = std::remove_if(neighbors.begin(), neighbors.end(), condition);
            neighbors.erase(it, neighbors.end());
        }
        
        void remove_duplicate_neighbors()
        {
            //defines a quasi-lexicographical ordering among set of neighbor objects
            const auto sort_predicate = [](const auto& n0, const auto& n1)
            {
                bool i1 = (n0.endpoint < n1.endpoint);
                auto lexi = 3 + 0*n0.edge;
                bool i0 = ctrs::dot_prod(lexi, 1+n0.edge) < ctrs::dot_prod(lexi, 1+n1.edge);
                return i1 || (!i1 && i0);
            };
            
            
            //defines an equivalence relation between two neighbor relationships
            const auto uniq_predicate = [](const auto& n0, const auto& n1)
            {
                return (n0.endpoint == n1.endpoint) && (n0.edge == n1.edge);
            };
            std::sort(neighbors.begin(), neighbors.end(), sort_predicate);
            auto it = std::unique(neighbors.begin(), neighbors.end(), uniq_predicate);
            neighbors.erase(it, neighbors.end());
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
                print(i.endpoint->terminal());
                print("----------------------");
            }
        }
        
        amr_node_t(){}
    };
}