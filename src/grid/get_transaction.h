#pragma once

#include "grid/transactions.h"
#include "amr/amr_node.h"

namespace spade::grid
{
    template <typename grid_t>
    auto get_transaction(const grid_t& grid, const std::size_t&, const neighbor_relation_t& relation)
    {
        grid_rect_copy_t output;
        auto lb_ini  = utils::tag[partition::global](relation.lb_ini);
        auto lb_term = utils::tag[partition::global](relation.lb_term);
        
        output.rank_send = grid.get_partition().get_rank(lb_ini);
        output.rank_recv = grid.get_partition().get_rank(lb_term);
        
        output.source.min(3) = grid.get_partition().to_local(lb_ini).value;
        output.source.max(3) = output.source.min(3)+1;
        output.dest.min(3)   = grid.get_partition().to_local(lb_term).value;
        output.dest.max(3)   = output.dest.min(3)+1;
        
        for (int d = 0; d < 3; ++d)
        {
            const int sign = relation.edge[d];
            const int nx   = grid.get_num_cells(d);
            const int ng   = grid.get_num_exchange(d);
            switch (sign)
            {
                case -1:
                {
                    output.source.min(d) = 0;
                    output.source.max(d) = ng;
                    output.dest.min(d)   = nx;
                    output.dest.max(d)   = nx+ng;
                    break;
                }
                case  0:
                {
                    output.source.min(d) = 0;
                    output.source.max(d) = nx;
                    output.dest.min(d)   = 0;
                    output.dest.max(d)   = nx;
                    break;
                }
                case  1:
                {
                    output.source.min(d) = nx-ng;
                    output.source.max(d) = nx;
                    output.dest.min(d)   = -ng;
                    output.dest.max(d)   = 0;
                    break;
                }
            }
        }
        return output;
    };
    
    template <typename grid_t>
    auto get_transaction(const grid_t& grid, const std::size_t& lb_ini, const amr::amr_neighbor_t<grid_t::dim()>& relation)
    {
        patch_fill_t<grid_t::dim()> output;
        neighbor_relation_t base_relation;
        base_relation.edge = relation.edge;
        base_relation.lb_ini  = lb_ini;
        base_relation.lb_term = relation.endpoint.get().tag;
        output.patches = get_transaction(grid, lb_ini, base_relation);
        
        auto& self  = grid.get_blocks().enumerated_nodes[lb_ini].get();
        auto& neigh = relation.endpoint.get();
        
        
        
        return output;
    };
}