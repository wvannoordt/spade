#pragma once

typedef double real_t;
typedef cvdf::ctrs::array<real_t, 3> v3d;
typedef cvdf::ctrs::array<real_t, 5> v5d;
typedef cvdf::ctrs::array<int,    3> v3i;
typedef cvdf::ctrs::array<int,    4> v4i;
typedef cvdf::ctrs::array<cvdf::grid::cell_t<int>, 4> v4c;
typedef cvdf::fluid_state::prim_t<real_t> prim_t;
typedef cvdf::fluid_state::cons_t<real_t> cons_t;

static void get_profile(const auto& q, std::vector<double>& y, std::vector<double>& u)
{
    std::vector<int> counts;
    const auto& grid  = q.get_grid();
    const auto& group = grid.group();
    auto rg   = grid.get_range(cvdf::grid::cell_centered);
    auto ymin = grid.get_bounds().min(1);
    int  ny   = grid.get_num_cells(1)*grid.get_num_blocks(1);
    y.resize(ny, 0.0);
    u.resize(ny, 0.0);
    counts.resize(ny, 0);
    for (auto i: rg)
    {
        const v4c  ijk(i[0], i[1], i[2], i[3]);
        const auto x  = grid.get_comp_coords(ijk);
        prim_t qloc;
        for (auto n: range(0,5))
        {
            qloc[n] = q(n, i[0], i[1], i[2], i[3]);
        }
        const auto xp = grid.get_coords(ijk);
        const auto dy = grid.get_dx(1);
        int idx  = round((x[1]-0.5*dy-ymin)/dy);
        y[idx] += xp[1];
        u[idx] += qloc.u();
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