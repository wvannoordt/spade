#include <chrono>
#include "spade.h"

typedef double real_t;

void diffuse(auto& rhs, const auto& q)
{
    const auto& grid = q.get_grid();
}

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
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
    real_t dt = 1e-3;
    
    auto calc_rhs = [&](auto& rhs, auto& q, const auto& t) -> void
    {
        rhs = 0.0;
        grid.exchange_array(q);
        diffuse(rhs, q);
    };
    
    spade::time_integration::rk2 time_int(phi, rhs, time0, dt, calc_rhs);
    
    spade::algs::fill_array(phi, ini);
    for (auto n: range(0, 10001))
    {
        print("nt = ", spade::utils::zfill(n, 8));
        time_int.advance();
        if (n%1000==0)
        {
            std::string fname = "phi";
            fname += spade::utils::zfill(n, 8);
            spade::io::output_vtk("output", fname, phi);
        }
    }
    return 0;
}
