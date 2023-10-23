#pragma once

#include <vector>

#include "core/ctrs.h"
#include "core/block_config.h"
#include "core/except.h"
#include "amr/amr_node.h"
#include "amr/amr_constraints.h"

namespace spade::amr
{    
    template <typename coord_val_t, typename array_designator_t>
    struct amr_blocks_t
    {
        using node_type       = amr::amr_node_t<array_designator_t::size()>;
        using handle_type     = node_handle_t<array_designator_t::size()>;
        using refine_type     = ctrs::array<bool, array_designator_t::size()>;
        using coord_val_type  = coord_val_t;
        
        ctrs::array<int, 3>                      num_blocks;
        bound_box_t<coord_val_t, 3>              bounds;
        std::vector<node_type>                   root_nodes;
        std::vector<bound_box_t<coord_val_t, 3>> block_boxes;
        std::vector<handle_type>                 all_nodes;
        std::vector<handle_type>                 enumerated_nodes;
        
        
        constexpr static std::size_t dim() { return array_designator_t::size(); }
        
        amr_blocks_t(const array_designator_t& num_blocks_in,
            const bound_box_t<coord_val_t, array_designator_t::size()>& bounds_in)
        {
            init(num_blocks_in, bounds_in);
        }
        
        void init(const auto& num_blocks_in, const auto& bounds_in)
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
                    if constexpr (dim() == 2) dijk[2] = 0;
                    bool all_zero = (dijk == 0*dijk);
                    if (!all_zero)
                    {
                        auto ijk_neigh = dijk + ijk;
                        for (auto k: range(0, dim()))
                        {
                            if (ijk_neigh[k]<0)              ijk_neigh[k] += num_blocks[k];
                            if (ijk_neigh[k]>=num_blocks[k]) ijk_neigh[k] -= num_blocks[k];
                        }
                        auto j = collapse_index(ijk_neigh, num_blocks);
                        
                        ctrs::array<int, dim()> dijk_d;
                        for (int d = 0; d < dim(); ++d) dijk_d[d] = dijk[d];
                        //node i has neighbor j
                        amr::amr_neighbor_t<dim()> neighbor_relationship{node_handle_t<dim()>(root_nodes, j), dijk_d};
                        root_nodes[i].neighbors.push_back(neighbor_relationship);
                    }
                }
                root_nodes[i].parent = handle_type::null();
                for (auto d: range(0, dim()))
                {
                    root_nodes[i].amr_position.min(d) = amr::amr_coord_t(ijk[d],   num_blocks[d], 0);
                    root_nodes[i].amr_position.max(d) = amr::amr_coord_t(ijk[d]+1, num_blocks[d], 0);
                }
                root_nodes[i].level = 0;
            }
            this->enumerate();
        }
        
        amr_blocks_t(const amr_blocks_t& rhs)
        {
            init(rhs.num_blocks, rhs.bounds);
            if (rhs.total_num_blocks() > rhs.num_blocks[0]*rhs.num_blocks[1]*rhs.num_blocks[2])
            {
                throw except::sp_exception("amr_blocks_t copy constructor not implemented for refined grid");
            }
        }
        
        auto& get_amr_node(const std::size_t lb_glob)
        {
            return enumerated_nodes[lb_glob].get();
        }
        
        const auto& get_amr_node(const std::size_t lb_glob) const
        {
            return enumerated_nodes[lb_glob].get();
        }
        
        //Note: these return the block bounding box in computational coordinates
        const auto& get_block_box (const std::size_t& lb_glob) const
        {
            return block_boxes[lb_glob];
        }
        
        std::size_t total_num_blocks() const { return enumerated_nodes.size(); }
        template <typename filter_t> std::size_t get_num_blocks(const filter_t& filter) const
        {
            std::size_t output = 0;
            for (const auto p: enumerated_nodes)
            {
                if (filter(*p)) ++output;
            }
            return output;
        }
        
        bound_box_t<bool, 3> is_domain_boundary(const std::size_t& i) const
        {
            return enumerated_nodes[i].get().is_domain_boundary();
        }
        
        // re-collects the global list of amr nodes and their bounding boxes,
        // must be called after every single refinement operation
        template <typename crit_t> void enumerate(const crit_t& criterion)
        {
            std::size_t block_count = 0;
            for (auto& n: root_nodes)
            {
                block_count += n.count_nodes();
            }
            all_nodes.clear();
            all_nodes.reserve(block_count);
            for (int i = 0; i < root_nodes.size(); ++i)
            {
                all_nodes.push_back(handle_type(root_nodes, i));
                root_nodes[i].collect_nodes(all_nodes);
            }
            
            
            enumerated_nodes.clear();
            enumerated_nodes.reserve(block_count);
            for (const auto& n: all_nodes)
            {
                if (criterion(n.get())) enumerated_nodes.push_back(n);
            }
            block_boxes.clear();
            block_boxes.resize(enumerated_nodes.size());
            ctrs::array<coord_val_t, dim()> dx;
            for (auto i: range(0, dim())) dx[i] = bounds.size(i)/num_blocks[i];
            for (auto i: range(0, block_boxes.size()))
            {
                const auto& amr_bbox  = enumerated_nodes[i].get().amr_position;
                auto& comp_bbox = block_boxes[i];
                comp_bbox.min(2) = 0.0;
                comp_bbox.max(2) = 1.0;
                for (auto d: range(0, dim()))
                {
                    comp_bbox.min(d) = amr_bbox.min(d).convert_to_coordinate(bounds.min(d), dx[d]);
                    comp_bbox.max(d) = amr_bbox.max(d).convert_to_coordinate(bounds.min(d), dx[d]);
                }
            }
            for (auto& p: all_nodes) p.get().tag = -1;
            for (int lb = 0; lb < enumerated_nodes.size(); ++lb) enumerated_nodes[lb].get().tag = lb;
        }
        
        void enumerate()
        {
            this->enumerate([](const auto& n){return n.terminal();});
        }
        
        template <typename interface_constraint_t>
        void refine_recurse(
            std::vector<handle_type>& lbs,
            const ctrs::array<bool, dim()>& is_periodic,
            const std::vector<typename node_type::amr_refine_t>& directions,
            const interface_constraint_t& constraint)
        {            
            std::vector<handle_type> new_children;
            int idx = 0;
            for (auto& lb: lbs)
            {
                auto& anode = lb.get();
                bool refd = anode.refine_node(directions[idx], is_periodic);
                if (refd)
                {
                    for (int i = 0; i < anode.subnodes.size(); ++i)
                    {
                        new_children.push_back(handle_type(anode.subnodes, i));
                    }
                }
                ++idx;
            }
            std::vector<handle_type>                      further_refines;
            std::vector<typename node_type::amr_refine_t> further_direcs;
            for (auto& child: new_children)
            {
                for (auto& neigh: child.get().neighbors)
                {
                    const auto idomain_boundary = child.get().is_domain_boundary();
                    auto new_refine = constraint(child.get(), neigh);
                    bool any = false;
                    for (auto b: new_refine) any = any | b;
                    bool reject_domain_boundary = false;
                    for (int i = 0; i < dim(); ++i)
                    {
                        reject_domain_boundary = reject_domain_boundary || ((neigh.edge[i] == -1) && idomain_boundary.min(i) && !is_periodic[i]);
                        reject_domain_boundary = reject_domain_boundary || ((neigh.edge[i] ==  1) && idomain_boundary.max(i) && !is_periodic[i]);
                    }
                    if (any && !reject_domain_boundary)
                    {
                        further_refines.push_back(neigh.endpoint);
                        further_direcs.push_back(new_refine);
                    }
                }
            }
            if (further_refines.size() > 0) this->refine_recurse(further_refines, is_periodic, further_direcs, constraint);
        }
        
        template <typename interface_constraint_t>
        void refine(
            std::vector<handle_type>& lbs,
            const ctrs::array<bool, dim()>& is_periodic,
            const std::vector<typename node_type::amr_refine_t>& directions,
            const interface_constraint_t& constraint)
        {
            this->refine_recurse(lbs, is_periodic, directions, constraint);
            this->enumerate();
        }
        
        template <typename interface_constraint_t>
        void refine(
            std::vector<handle_type>& lbs,
            const ctrs::array<bool, dim()>& is_periodic,
            const typename node_type::amr_refine_t& directions,
            const interface_constraint_t& constraint)
        {
            std::vector<typename node_type::amr_refine_t> drs(lbs.size(), directions);
            this->refine(lbs, is_periodic, drs, constraint);
        }
        
        template <typename interface_constraint_t>
        void refine(
            const std::vector<std::size_t>& lbs,
            const ctrs::array<bool, dim()>& is_periodic,
            const std::vector<typename node_type::amr_refine_t>& directions,
            const interface_constraint_t& constraint)
        {
            std::vector<handle_type> nds;
            for (auto lb: lbs)
            {
                nds.push_back(this->enumerated_nodes[lb]);
            }
            this->refine(nds, is_periodic, directions, constraint);
        }
        
        template <typename interface_constraint_t>
        void refine(std::size_t lb, const ctrs::array<bool, dim()>& is_periodic, const node_type::amr_refine_t directions, const interface_constraint_t& cc)
        {
            using v0_t = std::vector<std::size_t>;
            this->refine(v0_t{lb}, is_periodic, {directions}, cc);
        }
        
        template <typename interface_constraint_t>
        void refine(std::vector<std::size_t> lb, const ctrs::array<bool, dim()>& is_periodic, const node_type::amr_refine_t directions, const interface_constraint_t& cc)
        {
            using v0_t = std::vector<typename node_type::amr_refine_t>;
            v0_t v0;
            for (auto i:lb) v0.push_back(directions);
            this->refine(lb, is_periodic, v0, cc);
        }
        
        void refine(std::size_t lb, const ctrs::array<bool, dim()>& is_periodic, const node_type::amr_refine_t directions)
        {
            this->refine({lb}, is_periodic, {directions}, [](const auto&, const auto&){return typename node_type::amr_refine_t(false);});
        }
        
        void refine(const std::size_t& lb)
        {
            const ctrs::array<bool, dim()> is_per = false;
            const typename node_type::amr_refine_t directions = true;
            this->refine({lb}, is_per, directions);
        }
        
        const auto&       get_bounds()                                        const { return bounds; }
        const auto&       get_bounding_box(const std::size_t lb)              const { return block_boxes[lb]; }
        const coord_val_t get_size(const std::size_t i, const std::size_t lb) const { return block_boxes[lb].size(i); }
        const auto&       get_neighs(const std::size_t lb)                    const { return enumerated_nodes[lb].get().neighbors; }
        
        template <typename condition_t>
        std::vector<handle_type> select(const condition_t& condition)
        {
            std::vector<handle_type> output;
            for (auto& lb: enumerated_nodes)
            {
                if (condition(lb.get())) output.push_back(lb);
            }
            return output;
        }
    };
}