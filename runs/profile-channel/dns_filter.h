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
    
    static void extract_profile(profr_t& prof, const auto& q_o, const auto& q_i, const auto& callable)
    {
        std::vector<int> counts;
        const auto& grid  = q_o.get_grid();
        const auto& group = grid.group();
        if (group.isroot())
        {
            print("Extract", prof.name);
        }
        auto rg = grid.get_range(cvdf::grid::cell_centered);
        auto ymin = grid.get_bounds().min(1);
        int  ny   = grid.get_num_cells(1)*grid.get_num_blocks(1);
        
        counts.resize(ny, 0);
        for (auto i: rg)
        {
            const v4c  ijk(i[0], i[1], i[2], i[3]);
            const auto x  = grid.get_comp_coords(ijk);
            prim_t q_i_l, q_o_l;
            for (auto n: range(0,5))
            {
                q_i_l[n[0]] = q_i(n[0], i[0], i[1], i[2], i[3]);
                q_o_l[n[0]] = q_o(n[0], i[0], i[1], i[2], i[3]);
            }
            const auto xp = grid.get_coords(ijk);
            const auto dy = grid.get_dx(1);
            int idx  = round((x[1]-0.5*dy-ymin)/dy);
            prof.inst[idx] += callable(xp, q_o_l, q_i_l);
            counts[idx]++;
        }
        for (int ii = 0; ii < ny; ++ii)
        {
            prof.inst[ii] = group.sum(prof.inst[ii]);
            counts[ii]    = group.sum(counts[ii]);
        }
        for (int ii = 0; ii < ny; ++ii)
        {
            prof.inst[ii] /= counts[ii];
        }
    }
}
