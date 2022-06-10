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
    
    static void dns_filter(const cvdf::ctrs::array<int, 3>& filtsize, const auto& src, auto& dest)
    {
        const auto& grid = src.get_grid();
        const auto rg   = grid.get_range(cvdf::grid::cell_centered);
        v3i big;
        for (auto i: range(0,3)) big[i] = grid.get_num_cells(i)/filtsize[i];
        for (auto lb: range(0,grid.get_num_local_blocks()))
        {
            for (auto i: range(0,big[0])*range(0,big[1])*range(0,big[2]))
            {
                v3i ll;
                for (auto j: range(0, 3)) ll[j] = i[j]*filtsize[j];
                auto cell_range = range(ll[0],ll[0]+filtsize[0])*range(ll[1],ll[1]+filtsize[1])*range(ll[2],ll[2]+filtsize[2]);
                for (auto var: range(0, 5))
                {
                    real_t avg_val = 0.0;
                    for (auto c: cell_range) avg_val += src(var, c[0], c[1], c[2], lb);
                    avg_val /= cell_range.size();
                    for (auto c: cell_range) dest(var, c[0], c[1], c[2], lb) = avg_val;
                }
            }
        }
    }
    
    void noslip(const int& numcell, auto& prims)
    {
        const real_t t_wall = 0.1;
        const auto& grid = prims.get_grid();
        for (auto lb: range(0, grid.get_num_local_blocks()))
        {
            const auto& lb_glob = grid.get_partition().get_global_block(lb);
            int idc = 0;
            for (int dir = 2; dir <= 3; ++dir)
            {
                if (grid.is_domain_boundary(lb_glob, dir))
                {
                    const auto lb_idx = cvdf::ctrs::expand_index(lb_glob, grid.get_num_blocks());
                    const auto nvec_out = v3i(0,2*idc-1,0);
                    const cvdf::grid::cell_t<int> j = idc*(grid.get_num_cells(1)-1);
                    auto r1 = range(-grid.get_num_exchange(0), grid.get_num_cells(0) + grid.get_num_exchange(0));
                    auto r2 = range(-grid.get_num_exchange(2), grid.get_num_cells(2) + grid.get_num_exchange(2));
                    for (auto ii: r1*r2)
                    {
                        for (int nnn = 0; nnn < numcell; ++nnn)
                        {
                            v4c i_d(ii[0], j-(nnn+0)*nvec_out[1], ii[1], lb);
                            v4c i_g(ii[0], j+(nnn+1)*nvec_out[1], ii[1], lb);
                            prim_t q_d, q_g;
                            for (auto n: range(0,5)) q_d[n] = prims(n, i_d[0], i_d[1], i_d[2], i_d[3]);
                            const auto x_g = grid.get_comp_coords(i_g[0], i_g[1], i_g[2], i_g[3]);
                            const auto x_d = grid.get_comp_coords(i_d[0], i_d[1], i_d[2], i_d[3]);
                            const auto n_g = calc_normal_vector(grid.coord_sys(), x_g, i_g, 1);
                            const auto n_d = calc_normal_vector(grid.coord_sys(), x_d, i_d, 1);
                            q_g.p() =  q_d.p();
                            q_g.u() = -q_d.u();
                            q_g.v() = -q_d.v()*n_d[1]/n_g[1];
                            q_g.w() = -q_d.w();
                            q_g.T() =  t_wall;
                            for (auto n: range(0,5)) prims(n, i_g[0], i_g[1], i_g[2], i_g[3]) = q_g[n];
                        }
                    }
                }
                ++idc;
            }
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
        auto rg   = grid.get_range(cvdf::grid::cell_centered);
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
                q_i_l[n] = q_i(n, i[0], i[1], i[2], i[3]);
                q_o_l[n] = q_o(n, i[0], i[1], i[2], i[3]);
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
    
    static void extract_profile_f(profr_t& prof, const auto& q_o, const auto& q_i, const auto& callable)
    {
        std::vector<int> counts;
        const auto& grid  = q_o.get_grid();
        const auto& group = grid.group();
        if (group.isroot())
        {
            print("Extract", prof.name);
        }
        auto rg   = grid.get_range(cvdf::grid::cell_centered);
        auto ymin = grid.get_bounds().min(1);
        int  ny   = grid.get_num_cells(1)*grid.get_num_blocks(1);
        
        counts.resize(ny, 0);
        for (auto i: rg)
        {
            v4c  ijk_L(i[0], i[1], i[2], i[3]);
            v4c  ijk_R(i[0], i[1], i[2], i[3]);
            ijk_R[1] += 1;
            v5f  ijk_F(1, (int)ijk_L[0], (int)ijk_L[1], (int)ijk_L[2], (int)ijk_L[3]);
            const auto x  = grid.get_comp_coords(ijk_F);
            prim_t q_i_l_L, q_i_l_R, q_o_l_L, q_o_l_R;
            for (auto n: range(0,5))
            {
                q_i_l_L[n] = q_i(n, ijk_L[0], ijk_L[1], ijk_L[2], ijk_L[3]);
                q_o_l_L[n] = q_o(n, ijk_L[0], ijk_L[1], ijk_L[2], ijk_L[3]);
                q_i_l_R[n] = q_i(n, ijk_R[0], ijk_R[1], ijk_R[2], ijk_R[3]);
                q_o_l_R[n] = q_o(n, ijk_R[0], ijk_R[1], ijk_R[2], ijk_R[3]);
            }
            const auto xp = grid.get_coords(ijk_F);
            const auto dy = grid.get_dx(1);
            int idx  = round((xp[1]-ymin)/dy);
            if (idx < 0 || idx >= ny)
            {
                print(idx, ijk_F, xp, ymin, dy);
                group.pause();
            }
            prof.inst[idx] += callable(xp, q_o_l_L, q_i_l_L, q_o_l_R, q_i_l_R);
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
