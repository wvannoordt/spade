#pragma once

#include "core/except.h"
#include "core/ctrs.h"
#include "grid/transactions.h"
#include "amr/amr_node.h"

namespace spade::grid
{
    template <typename nexchg_t, typename grid_t>
    auto get_transaction(
        const nexchg_t& num_exchgs,
        const grid_t& src_grid,
        const grid_t& dst_grid,
        const std::size_t&,
        const neighbor_relation_t& relation)
    {
        grid_rect_copy_t output;
        auto lb_ini  = utils::tag[partition::global](relation.lb_ini);
        auto lb_term = utils::tag[partition::global](relation.lb_term);
        
        output.rank_send = src_grid.get_partition().get_rank(lb_ini);
        output.rank_recv = dst_grid.get_partition().get_rank(lb_term);
        
        output.glob_source_blk = lb_ini.value;
        output.glob_dest_blk   = lb_term.value;

        output.source.min(3) = src_grid.get_partition().to_local(lb_ini);
        output.source.max(3) = output.source.min(3) + 1;
        output.dest.min(3)   = dst_grid.get_partition().to_local(lb_term);
        output.dest.max(3)   = output.dest.min(3) + 1;
        
        // Here we need to classify this transaction
        int dsum = 0;
        for (const auto& ii: relation.edge) dsum += utils::abs(ii);
        ctrs::array<transaction_tag_t, 4> tags{null_transaction, xface_transaction, xedge_transaction, crnr_transaction};
        int tag_id = int(tags[dsum]);
        if (dsum == 1) // face
        {
            for (int d = 0; d < relation.edge.size(); ++d)
            {
                if (relation.edge[d] != 0) tag_id -= d;
            }
        }
        if (dsum == 2) // edge
        {
            for (int d = 0; d < relation.edge.size(); ++d)
            {
                if (relation.edge[d] == 0) tag_id -= d;
            }
        }
        output.tag = transaction_tag_t(tag_id);
        // Transaction tag has been computed.
        
        for (int d = 0; d < 3; ++d)
        {
            const int sign = relation.edge[d];
            const int nx   = src_grid.get_num_cells(d);
            const int ng   = num_exchgs[d];
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
        if constexpr (grid_t::dim() == 2)
        {
            output.source.min(2) = 0;
            output.source.max(2) = 1;
            output.dest.min(2)   = 0;
            output.dest.max(2)   = 1;
        }
        
        return output;
    };
    
    template <typename nexchg_t, typename grid_t>
    auto get_transaction(
        const nexchg_t& num_exchgs,
        const grid_t& src_grid,
        const grid_t& dst_grid,
        const std::size_t& lb_ini,
        const amr::amr_neighbor_t<grid_t::dim()>& relation)
    {
        patch_fill_t<grid_t::dim()> output;
        neighbor_relation_t base_relation;
        base_relation.edge = -1*relation.edge;
        if constexpr (grid_t::dim() == 2) base_relation.edge[2] = 0;
        base_relation.lb_ini  = relation.endpoint.get().tag;
        base_relation.lb_term = lb_ini;
        
        
        //NOTE: here, we ASK the neighbor for data instead of tell it about the data we send.
        auto ptch = get_transaction(num_exchgs, src_grid, dst_grid, lb_ini, base_relation);
        output.patches = ptch;
        
        auto& self  = relation.endpoint.get();
        auto& neigh = dst_grid.get_blocks().get_amr_node(lb_ini);
        
        int tag = int(ptch.tag);
        
        for (int d = 0; d < grid_t::dim(); ++d)
        {
            int level_diff = self.level[d] - neigh.level[d];
            int edgediff   = base_relation.edge[d];
            if (edgediff == 0 && level_diff > 0)
            {
                tag += 3;
            }
        }
        
        output.tag = transaction_tag_t(tag);
        
        output.i_coeff = 0;
        output.i_incr  = 0;
        for (int d = 0; d < grid_t::dim(); ++d)
        {
            int level_diff = self.level[d] - neigh.level[d];
            
            //Need to modify:
            
            // [ ] output.patches.source
            // [ ] output.patches.dest
            // [x] output.i_skip
            // [x] output.num_oslot
            // [x] output.delta_i
            
            switch(level_diff)
            {
                case -1: // neighbor is finer, I am coarser in this direction
                {
                    output.i_coeff[d] = -1;
                    // Neighbor is finer, so I only need half the range in this
                    // direction
                    switch(base_relation.edge[d])
                    {
                        case -1:
                        {
                            output.patches.source.max(d) = output.patches.source.min(d) + output.patches.source.size(d)/2;
                            break;
                        }
                        case  0:
                        {
                            output.patches.source.max(d) = output.patches.source.min(d) + output.patches.source.size(d)/2;
                            if (neigh.amr_position.min(d)>self.amr_position.min(d))
                            {
                                const auto sz = output.patches.source.size(d);
                                output.patches.source.min(d) += sz;
                                output.patches.source.max(d) += sz;
                            }
                            break;
                        }
                        case  1:
                        {
                            output.patches.source.min(d) = output.patches.source.max(d) - output.patches.source.size(d)/2;
                            break;
                        }
                    }
                    break;
                }
                case  0: // same level (do nothing)
                {
                    break;
                }
                case  1: // neighbor is coarser, I am finer in this direction
                {
                    output.i_incr[d] = 1;
                    output.i_coeff[d] = 1;
                    switch(base_relation.edge[d])
                    {
                        // neighbor is coarser, so send twice the data
                        case -1:
                        {
                            output.patches.source.max(d) = output.patches.source.min(d) + 2*output.patches.source.size(d);
                            break;
                        }
                        case  0:
                        {
                            //except here, where we only fill half the data!
                            output.patches.dest.max(d) = output.patches.dest.min(d) + output.patches.dest.size(d)/2;
                            if (self.amr_position.min(d) > neigh.amr_position.min(d))
                            {
                                const auto sz = output.patches.dest.size(d);
                                output.patches.dest.min(d) += sz;
                                output.patches.dest.max(d) += sz;
                            }
                            break;
                        }
                        case  1:
                        {
                            output.patches.source.min(d) = output.patches.source.max(d) - 2*output.patches.source.size(d);
                            break;
                        }
                    }
                    
                    break;
                }
                default:
                {
                    std::stringstream message;
                    message << "Found illegal grid configuration.\n";
                    message << "Self.level:  " << self.level  << "\n";
                    message << "Other.level: " << neigh.level << "\n";
                    // message << ""
                    
                    const auto bbx  = dst_grid.get_bounding_box(utils::tag[partition::global](lb_ini));
                    const auto bbx2 = src_grid.get_bounding_box(utils::tag[partition::global](base_relation.lb_ini));
                    using bbx_t = typename utils::remove_all<decltype(bbx)>::type;
                    using r_t = typename bbx_t::value_type;
                    constexpr std::size_t dim = bbx.size();
                    std::vector<bbx_t> bbxs{bbx, bbx2};
                    throw except::bbxs_exception<r_t, dim>(message.str(), bbxs);
                }
            }
        }
        
        // if (global::debug)
        // {
        //     const auto& group = src_grid.group();
        //     int lb_loc = output.patches.dest.min(3);
        //     int lb_glob = src_grid.get_partition().local_block_to_global_block[lb_loc];
        //     group.exclusive([&](){
        //     if (lb_glob == 6859 && relation.edge[0] == 0 && relation.edge[1] == 0 && relation.edge[2] < 0)
        //     {
        //         print("=========================");
        //         print("lb loc/glob", lb_loc, lb_glob);
        //         print("rank", group.rank());
        //         print("reducible", output.reducible());
        //         print("source", output.patches.source, src_grid.get_partition().local_block_to_global_block[output.patches.source.min(3)]);
        //         print("dest", output.patches.dest, src_grid.get_partition().local_block_to_global_block[output.patches.dest.min(3)]);
        //         print("sender", output.sender());
        //         print("receiver", output.receiver());
        //         print("lb_ini (G)", lb_ini);
        //         print("lb_term (G)", relation.endpoint->tag);
        //         print("=========================");
        //     }
        //     });
        // }
        
        return output;
    };
}