#pragma once

#include <vector>
#include <tuple>
#include <algorithm>

#include "core/utils.h"
#include "core/ctrs.h"
#include "core/bounding_box.h"
#include "amr/amr_coord.h"

namespace spade::amr
{    
    template <const std::size_t grid_dim, typename tag_t = int> struct amr_node_t;
    template <const std::size_t grid_dim> using node_handle_t = utils::vector_location<amr_node_t<grid_dim>>;
    template <const std::size_t grid_dim> struct amr_neighbor_t
    {
        //                            edge
        //         amr_node          ------>     endpoint
        node_handle_t<grid_dim>  endpoint;
        ctrs::array<int, grid_dim> edge;
    };
    template <const std::size_t grid_dim, typename tag_t> struct amr_node_t
    {
        using amr_refine_t = ctrs::array<bool,grid_dim>;
        using handle_type  = node_handle_t<grid_dim>;
        
        bound_box_t<amr_coord_t, grid_dim>    amr_position;
        std::vector<amr_node_t>               subnodes;
        std::vector<amr_neighbor_t<grid_dim>> neighbors;
        ctrs::array<int, grid_dim>            level;
        handle_type                           parent = handle_type::null();
        tag_t                                 tag; //this holds the blocks value!!!!!!!
        
        amr_node_t(const amr_node_t&) = default;
        
        static constexpr std::size_t dim() { return grid_dim; }
        
        bool terminal() const { return subnodes.size()==0; }
        
        std::size_t count_nodes() const
        {
            std::size_t output = 1;
            for (auto& n: subnodes) output += n.count_nodes();
            return output;
        }
        
        bool refine_node(const amr_refine_t& ref_type, const ctrs::array<bool, dim()> is_periodic = false)
        {
            std::size_t num_sub_nodes = 1;
            for (auto i: range(0, dim())) num_sub_nodes *= ref_type[i]?2:1;
            
            //exception conditions
            if (num_sub_nodes == 1) return false;
            if (!this->terminal())  return false;
            subnodes.resize(num_sub_nodes);
            auto rg = range(0,ref_type[0]?2:1)*range(0,ref_type[1]?2:1)*range(0,(ref_type[dim()-1]&&(dim()==3))?2:1);
            int ct = 0;
            for (auto i: rg)
            {
                auto& child = subnodes[ct];
                child.parent = handle_type(subnodes, ct);
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
                handle_type child_handle(subnodes, i);
                
                //First add neighbor relationships for siblings of current child
                for (int j = 0; j < subnodes.size(); ++j)
                {
                    // if (i != j)
                    // {
                        handle_type sibling(subnodes, j);
                        child.create_neighbor(sibling);
                    // }
                }
                
                //Now add parent neighbors
                for (auto& neigh: neighbors)
                {
                    child.create_neighbor(neigh.endpoint);
                    neigh.endpoint.get().create_neighbor(child_handle);
                }
            }
            
            //Each of my neighbors track their own neighbors, so need to update those.
            //Now that I am no longer a terminal node, remove myself from my neighbors' list of neighbors
            for (auto& neigh: neighbors) neigh.endpoint.get().remove_neighbors([&](const auto& ngh) { return !ngh.endpoint.get().terminal(); });

            for (auto& child: subnodes) child.remove_duplicate_neighbors();
            for (auto& neigh: neighbors) neigh.endpoint.get().remove_duplicate_neighbors();
            
            return true;
        }
        
        template <typename condition_t>
        void create_neighbor(handle_type& candidate, const condition_t& condition, const bool periodic_neighs)
        {
            if (!candidate.get().terminal()) return;
            ctrs::array<ctrs::array<bool, 3>, 3> table;
            table.fill(true);
            
            for (int d = 0; d < dim(); ++d)
            {
                auto& in_region = table[d];
                const auto& n_min = this->amr_position.min(d);
                const auto& n_max = this->amr_position.max(d);
                
                const auto& m_min = candidate.get().amr_position.min(d);
                const auto& m_max = candidate.get().amr_position.max(d);
                
                const auto uboundary_wrap = [periodic_neighs](const auto& test, const auto& x0, const auto& x1)
                {
                    auto y0 = x0;
                    auto y1 = x1;
                    if ((test.partition == test.num_partitions) && periodic_neighs)
                    {
                        y1.partition += y1.num_partitions;
                        y0.partition += y0.num_partitions;
                    }
                    return std::make_tuple(y0, y1);
                };
                
                const auto lboundary_wrap = [periodic_neighs](const auto& test, const auto& x0, const auto& x1)
                {
                    auto y0 = x0;
                    auto y1 = x1;
                    if (((test.partition == 0) && (test.bits == 0)) && periodic_neighs)
                    {
                        y1.partition -= y1.num_partitions;
                        y0.partition -= y0.num_partitions;
                    }
                    return std::make_tuple(y0, y1);
                };
                
                //Wrap around boundaries to detect periodic neighbors
                const auto [m_min_m, m_max_m] = lboundary_wrap(n_min, m_min, m_max);
                const auto [m_min_p, m_max_p] = uboundary_wrap(n_max, m_min, m_max);
                
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
                        bool allzero = edge == 0;
                        if (valid && !allzero)
                        {
                            amr_neighbor_t<dim()> neighbor{candidate, 0};
                            ctrs::copy_array(edge, neighbor.edge);
                            if (condition(neighbor)) neighbors.push_back(neighbor);
                        }
                    }
                }
            }
        }
        
        void create_neighbor(handle_type& candidate)
        {
            this->create_neighbor(candidate, [](const auto&){return true;}, true);
        }
        
        void create_neighbor(handle_type& candidate, const bool periodic_neighs)
        {
            this->create_neighbor(candidate, [](const auto&){return true;}, periodic_neighs);
        }
        
        template <typename condition_t>
        void remove_neighbors(const condition_t& condition)
        {
            auto it = std::remove_if(neighbors.begin(), neighbors.end(), condition);
            neighbors.erase(it, neighbors.end());
        }
        
        void remove_duplicate_neighbors()
        {
            //If you are making changes to this function, then BE VERY CAREFUL!!
            
            //defines a quasi-lexicographical ordering among set of neighbor objects
            const auto sort_comp = [](const auto& n0, const auto& n1)
            {
                // if (&n0.endpoint.get() == &n1.endpoint.get())
                // {
                //     auto lexi = 0*n0.edge;
                //     lexi[0] = 1;
                //     for (int i = 1; i < lexi.size(); ++i) lexi[i] = 3*lexi[i-1];
                //     return ctrs::dot_prod(lexi, 1+n0.edge) < ctrs::dot_prod(lexi, 1+n1.edge);
                // }
                // return &n0.endpoint.get() < &n1.endpoint.get();
                const auto& node0 = *n0.endpoint;
                const auto& node1 = *n1.endpoint;
                
                const auto& edge0 = n0.edge;
                const auto& edge1 = n1.edge;
                
                for (int d = 0; d < edge0.size(); ++d)
                {
                    if (edge0[d] != edge1[d]) return edge0[d] < edge1[d];
                }
                
                // Edge vectors are equal here
                const auto& amrbbx0 = node0.amr_position;
                const auto& amrbbx1 = node1.amr_position;
                for (int d = 0; d < amrbbx0.size(); ++d)
                {
                    if (amrbbx0.min(d) != amrbbx1.min(d)) return amrbbx0.min(d) < amrbbx1.min(d);
                    if (amrbbx0.max(d) != amrbbx1.max(d)) return amrbbx0.max(d) < amrbbx1.max(d);
                }
                
                //Everything is equal here...
                return amrbbx0.max(0) < amrbbx1.max(0);
            };
            
            
            //defines an equivalence relation between two neighbor relationships
            const auto uniq_predicate = [&](const auto& n0, const auto& n1)
            {
                return !sort_comp(n0, n1) && !sort_comp(n1, n0);
            };

            std::stable_sort(neighbors.begin(), neighbors.end(), sort_comp);
            auto it = std::unique(neighbors.begin(), neighbors.end(), uniq_predicate);
            neighbors.erase(it, neighbors.end());
        }
        
        void collect_nodes(std::vector<handle_type>& node_list)
        {
            for (int i = 0; i < subnodes.size(); ++i)
            {
                node_list.push_back(handle_type(subnodes, i));
            }
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
            out.bnds = false;
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
                print(i.edge, "@", &i.endpoint.get());
                print(i.endpoint.get().ubox());
                print(i.endpoint.get().terminal());
                print("----------------------");
            }
        }
        
        amr_node_t(){neighbors.reserve(100);}
    };
}