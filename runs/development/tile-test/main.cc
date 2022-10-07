#include <chrono>
#include "spade.h"

typedef double real_t;

void diffuse(auto& rhs, const auto& q)
{
    const auto& grid = q.get_grid();
    int nlb = grid.get_num_local_blocks();
    int nk = grid.get_num_cells(2);
    int nj = grid.get_num_cells(1);
    int ni = grid.get_num_cells(0);
    real_t rhsloc = 0.0;
    for (int lb = 0; lb < nlb; ++lb)
    {
        auto dx = grid.get_block_box(lb).size(0)/ni;
        auto dy = grid.get_block_box(lb).size(1)/nj;
        auto dz = grid.get_block_box(lb).size(2)/nk;
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    rhsloc  = (q(i+1,   j,   k, lb) - 2.0*q(i, j, k, lb) + q(i-1,   j,   k, lb))/(dx*dx);
                    rhsloc += (q(i,   j+1,   k, lb) - 2.0*q(i, j, k, lb) + q(i,   j-1,   k, lb))/(dy*dy);
                    rhsloc += (q(i,   j,   k+1, lb) - 2.0*q(i, j, k, lb) + q(i,   j,   k-1, lb))/(dz*dz);
                    rhs(i, j, k, lb) += rhsloc;
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3> cells_in_block(160, 160, 160);
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
    spade::utils::mtimer_t tmr;
    
    auto calc_rhs = [&](auto& rhs, auto& q, const auto& t) -> void
    {
        rhs = 0.0;
        grid.exchange_array(q);
        diffuse(rhs, q);
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
