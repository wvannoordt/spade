#pragma once

#include "ibm/boundary_info.h"

namespace spade::ibm
{
    template <typename ghosts_t, typename grid_t, typename geom_t, typename sample_dist_t>
    auto compute_ghost_sample_points(const ghosts_t& ghosts, const grid_t& grid, const geom_t& geom, const sample_dist_t& sfunc)
    {
        using pnt_t = typename ghosts_t::pnt_t;
        spade::device::shared_vector<pnt_t> ips;
        std::size_t idx = 0;
        for (int dir = 0; dir < ghosts.aligned.size(); ++dir)
        {
            std::size_t count = ghosts.aligned[dir].indices.size();
            for (std::size_t idx = 0; idx < count; ++idx)
            {
                for (int layer = 0; layer < ghosts.aligned[dir].num_layer(); ++layer)
                {
                    auto nv = ghosts.aligned[dir].closest_normals[idx][layer];
                    auto xc = ghosts.aligned[dir].closest_points[idx][layer];
                    // if (!ghosts.aligned[dir].can_fill[idx][layer])
                    // {
                    //     nv = ghosts.aligned[dir].closest_normals[idx][0];
                    //     xc = ghosts.aligned[dir].closest_points[idx][0];
                    // }
                    const auto sampldist = sfunc;
                    pnt_t ip = xc;
                    ip += sampldist*nv;
                    
                    if (geom.is_interior(ip))
                    {
                        auto irreg_idx = ghosts.aligned[dir].indices[idx][0];
                        const auto gp_sign = ghosts.aligned[dir].signs[idx];
                        irreg_idx.i(dir) += gp_sign;
                        ip = grid.get_coords(irreg_idx);
                    }
                    
                    ips.push_back(ip);
                }
            }
        }
        std::size_t diag_count = ghosts.diags.indices.size();
        for (std::size_t idx = 0; idx < diag_count; ++idx)
        {
            auto nv = ghosts.diags.closest_normals[idx];
            auto xc = ghosts.diags.closest_points[idx];
            const auto sampldist = sfunc;
            pnt_t ip = xc;
            ip += sampldist*nv;
            ips.push_back(ip);
        }
        ips.transfer();
        return ips;
    }
}