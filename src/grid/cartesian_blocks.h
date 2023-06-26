#pragma once

#include <vector>

#include "core/ctrs.h"
#include "core/bounding_box.h"

namespace spade::grid
{
    // this is a globally-owned arrangement of blocks.
    //
    // specifies the location in computational space
    // of each block
    //
    // Specifies the neighbour relationships of each block
    
    struct neighbor_relation_t
    {
        std::size_t lb_ini, lb_term;
        ctrs::array<int, 3> edge;
    };
    
    template <typename array_designator_t, typename float_t>
    struct cartesian_blocks_t
    {
        
        static constexpr std::size_t conf_dim = array_designator_t::size();
        static constexpr std::size_t dim() { return conf_dim; }
        
        std::size_t                           total_blocks;
        ctrs::array<int, 3>                   num_blocks;
        bound_box_t<float_t,  3>              bounds;
        std::vector<bound_box_t<float_t, 3>>  block_boxes;
        std::vector<bound_box_t<bool, 3>>     block_is_domain_boundary;
        ctrs::array<float_t, 3> bsize;
        std::vector<ctrs::array<neighbor_relation_t, static_math::pow<3, dim()>::value - 1>> neighbors;
        
        cartesian_blocks_t(const array_designator_t& num_blocks_in, const bound_box_t<float_t, array_designator_t::size()>& bounds_in)
        {
            //There is no notion of local and global here.
            ctrs::copy_array(num_blocks_in, num_blocks, 1);
            bounds.min(2) = 0.0;
            bounds.max(2) = 1.0;
            for (int i = 0; i < dim(); ++i)
            {
                bounds.min(i) = bounds_in.min(i);
                bounds.max(i) = bounds_in.max(i);
            }
            for (int i = 0; i < bsize.size(); ++i) bsize[i] = bounds.size(i)/num_blocks[i];
            total_blocks = num_blocks[0]*num_blocks[1]*num_blocks[2];
            block_is_domain_boundary.resize(total_blocks);
            block_boxes.resize(total_blocks);
            neighbors.resize(total_blocks);
            auto rg = range(0, num_blocks[0])*range(0, num_blocks[1])*range(0, num_blocks[2]);
            int lb = 0;
            for (auto blk: rg)
            {
                algs::static_for<0, 3>([&](const auto& ii)
                {
                    const int iv = ii.value;
                    const auto idx = udci::idx_const_t<iv>();
                    if (blk[idx] == 0)                block_is_domain_boundary[lb].min(iv) = true;
                    if (blk[idx] == num_blocks[iv]-1) block_is_domain_boundary[lb].max(iv) = true;
                });
                
                auto& bnd = block_boxes[lb];
                for (int i = 0; i < num_blocks.size(); ++i)
                {
                    bnd.min(i) = bounds.min(i) + blk[i]*bsize[i];
                    bnd.max(i) = bnd.min(i)    + bsize[i];
                }
                
                auto& neighs = neighbors[lb];
                int dzmin = -1, dzmax = 2;
                if constexpr (dim() == 2) { dzmin = 0; dzmax = 1; } 
                auto lb_rg = range(-1,2)*range(-1,2)*range(dzmin,dzmax);
                int neigh_idx = 0;
                for (auto edge: lb_rg)
                {
                    if (!((edge[0_c] == 0) && (edge[1_c] == 0) && (edge[2_c] == 0)))
                    {
                        
                        ctrs::array<int, 3> lb_neigh;
                        ctrs::array<int, 3> lb_pos(blk[0_c], blk[1_c], blk[2_c]);
                        algs::static_for<0, 3>([&](const auto& ii)
                        {
                            const int iv = ii.value;
                            const auto idx = udci::idx_const_t<iv>();
                            lb_neigh[iv] = blk[idx] + edge[idx];
                            if (lb_neigh[iv] <  0)              lb_neigh[iv] += num_blocks[iv];
                            if (lb_neigh[iv] >= num_blocks[iv]) lb_neigh[iv] -= num_blocks[iv];
                        });
                        
                        neighbor_relation_t new_data;
                        new_data.lb_ini  = lb;
                        new_data.lb_term = ctrs::collapse_index(lb_neigh, num_blocks);
                        new_data.edge[0] = edge[0_c];
                        new_data.edge[1] = edge[1_c];
                        new_data.edge[2] = edge[2_c];
                        neighbors[lb][neigh_idx] = new_data;
                        
                        ++neigh_idx;
                    }
                }
                ++lb;
            }
        }
        
        const auto& get_neighs(const std::size_t& i)             const { return neighbors[i]; }
        const auto& get_bounding_box(const std::size_t& i)       const { return block_boxes[i]; }
        const auto& is_domain_boundary(const std::size_t& i)     const { return block_is_domain_boundary[i]; }
        std::size_t total_num_blocks()                           const { return total_blocks; }
        const auto get_size(const std::size_t& i)                const { return bsize; }
        const auto get_size(const int dir, const std::size_t& i) const { return bsize[dir]; }
        const auto& get_bounds()                                 const { return bounds; }
    };
}