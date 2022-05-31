#pragma once

#include "cvdf.h"
#include "local_types.h"

namespace postprocessing
{
    static void copy_field(const auto& src, auto& dest)
    {
        const auto rg = 
            range(0, src.get_minor_dims().total_size())*
            src.get_grid().get_range(cvdf::grid::cell_centered)*
            range(0, src.get_major_dims().total_size());
        for (auto i: rg)
        {
            dest.unwrap_idx(i[0], i[1], i[2], i[3], i[4], i[5]) = src.unwrap_idx(i[0], i[1], i[2], i[3], i[4], i[5]);
        }
    }
    
    static void dns_filter(const auto& src, auto& dest)
    {
        const int filtnum = 8;
        const auto& grid = src.get_grid();
        const auto rg   = grid.get_range(cvdf::grid::cell_centered);
        const auto rg_f = range(-filtnum, filtnum+1)*range(-filtnum, filtnum+1)*range(-filtnum, filtnum+1);
        for (auto i: range(0,5)*rg)
        {
            dest(i[0], i[1], i[2], i[3], i[4]) = 0.0;
            for (auto di: rg_f)
            {
                dest(i[0], i[1], i[2], i[3], i[4]) += src(i[0], i[1]+di[0], i[2]+di[1], i[3]+di[2], i[4]);
            }
            dest(i[0], i[1], i[2], i[3], i[4]) /= (real_t)(rg_f.size());
        }
    }
    
    static void extract_vel_profile(const auto& q, std::vector<real_t>& y, std::vector<real_t>& u)
    {
        std::vector<int> counts;
        const auto& grid  = q.get_grid();
        const auto& group = grid.group();
        auto rg = grid.get_range(cvdf::grid::cell_centered);
        auto ymin = grid.get_bounds().min(1);
        int  ny   = grid.get_num_cells(1)*grid.get_num_blocks(1);
        
        counts.resize(ny, 0);
        y.resize(ny, 0.0);
        u.resize(ny, 0.0);
        for (auto i: rg)
        {
            const v4c  ijk(i[0], i[1], i[2], i[3]);
            const auto x  = grid.get_comp_coords(ijk);
            const auto xp = grid.get_coords(ijk);
            const auto dy = grid.get_dx(1);
            int idx  = round((x[1]-0.5*dy-ymin)/dy);
            y[idx] += xp[1];
            u[idx] += q(2, i[0], i[1], i[2], i[3]);
            counts[idx]++;
        }
        for (int ii = 0; ii < ny; ++ii)
        {
            y[ii]      = group.sum(y[ii]);
            u[ii]      = group.sum(u[ii]);
            counts[ii] = group.sum(counts[ii]);
        }
        for (int ii = 0; ii < ny; ++ii)
        {
            y[ii] /= counts[ii];
            u[ii] /= counts[ii];
        }
    }
}
