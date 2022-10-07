#include <chrono>
#include "spade.h"

typedef double real_t;

template <typename fcn_t> void ijk_loop(int ni, int nj, int nk, const fcn_t& fcn)
{
    for (int k = 0; k < nk; ++k)
    {
        for (int j = 0; j < nj; ++j)
        {
            for (int i = 0; i < ni; ++i)
            {
                fcn(i, j, k);
            }
        }
    }
}

void diffuse(auto& rhs, const auto& q)
{
    const auto& grid = q.get_grid();
    int nlb = grid.get_num_local_blocks();
    int nk = grid.get_num_cells(2);
    int nj = grid.get_num_cells(1);
    int ni = grid.get_num_cells(0);
    const real_t a0 = -1.0/12.0;
    const real_t a1 = 4.0/3.0;
    const real_t a2 = -5.0/2.0;
    const real_t a3 = 4.0/3.0;
    const real_t a4 = -1.0/12.0;
    for (int lb = 0; lb < nlb; ++lb)
    {
        auto dx = grid.get_block_box(lb).size(0)/ni;
        auto dy = grid.get_block_box(lb).size(1)/nj;
        auto dz = grid.get_block_box(lb).size(2)/nk;
        auto fcn = [&](int i, int j, int k) -> void
        {
        // for (int k = 0; k < nk; ++k)
        // for (auto k: range(0,nk))
        // {
            // for (int j = 0; j < nj; ++j)
            // for (auto j: range(0,nj))
            // {
                // for (int i = 0; i < ni; ++i)
                // for (auto i: range(0,ni))
                // {
        // for (auto idx: range(0,ni)*range(0,nj)*range(0,nk))
        // {
            // int i = idx[0];
            // int j = idx[1];
            // int k = idx[2];
                    int i = idx[0];
                    int j = idx[1];
                    int k = idx[2];
                    rhs(i, j, k, lb) += (a0*q(i+2,   j,   k, lb) + a1*q(i+1,   j,   k, lb) + a2*q(i, j, k, lb) + a3*q(i-1,   j,   k, lb) + a4*q(i-2,   j,   k, lb))/(dx*dx);
                    rhs(i, j, k, lb) += (a0*q(i,   j+2,   k, lb) + a1*q(i,   j+1,   k, lb) + a2*q(i, j, k, lb) + a3*q(i,   j-1,   k, lb) + a4*q(i,   j-2,   k, lb))/(dy*dy);
                    rhs(i, j, k, lb) += (a0*q(i,   j,   k+2, lb) + a1*q(i,   j,   k+1, lb) + a2*q(i, j, k, lb) + a3*q(i,   j,   k-1, lb) + a4*q(i,   j,   k-2, lb))/(dz*dz);
        // }
                // }
            // }
        // }
        };
        ctrs::array<int, 3> nijk(ni, nj, nk);
        cell_loop(nijk, fcn);
    }
}

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(20, 20, 20);
    spade::ctrs::array<int, 3> cells_in_block(16, 16, 16);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    
    bounds.min(0) = 0.0;
    bounds.max(0) = 1.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 1.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 1.0;
    
    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    real_t fill = 0.0;
    
    spade::grid::cell_array phi(grid, fill);
    spade::grid::cell_array rhs(grid, fill);
    
    spade::ctrs::array<real_t, 3> x0;
    for (auto i: range(0,3)) x0[i] = 0.5*(bounds.min(i) + bounds.max(i));
    auto ini = [&](const spade::ctrs::array<real_t, 3> x) -> real_t
    {
        spade::ctrs::array<real_t, 3> dx = x - x0;
        real_t r = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        return std::exp(-15.0*r*r);
    };
    
    real_t time0 = 0.0;
    real_t dt = 1e-5;
    spade::utils::mtimer_t tmr(3);
    
    auto calc_rhs = [&](auto& rhs, auto& q, const auto& t) -> void
    {
        tmr.start(0);
        rhs = 0.0;
        tmr.stop(0);
        
        tmr.start(1);
        grid.exchange_array(q);
        tmr.stop(1);
        
        tmr.start(2);
        diffuse(rhs, q);
        tmr.stop(2);
        
        print(tmr);
    };
    
    spade::time_integration::rk2 time_int(phi, rhs, time0, dt, calc_rhs);
    
    std::string frmt = "RHS: {}/{}    Exchange: {}/{}  Adv: {}/{}";
    spade::algs::fill_array(phi, ini);
    for (auto n: range(0, 10001))
    {
        time_int.advance();
        print("nt = ", spade::utils::zfill(n, 8));
        if (n%1000==0)
        {
            std::string fname = "phi";
            fname += spade::utils::zfill(n, 8);
            spade::io::output_vtk("output", fname, phi);
        }
    }
    return 0;
}
