#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<real_t, 3> v3d;
typedef spade::ctrs::array<int,    3> v3i;
typedef spade::ctrs::array<int,    4> v4i;
typedef spade::ctrs::array<spade::grid::cell_t<int>, 4> v4c;
typedef spade::fluid_state::prim_t<real_t> prim_t;
typedef spade::fluid_state::cons_t<real_t> cons_t;

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

    prim_t fill = 0.0;
    spade::grid::grid_array prim (grid, fill);
    
    auto ini = [&](const spade::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        prim_t output;
        output.p() = 0.5*x[0];
        output.T() = 0.5*x[1];
        output.u() = 0.5*sin(x[0]);
        output.v() = 0.9*cos(x[1]);
        output.w() = 1.2*cos(x[2]);
        return output;
    };
    
    spade::algs::fill_array(prim, ini);
    spade::io::output_vtk("output", "ini", grid, prim);
    spade::io::binary_write("test.bin", prim);
    prim = 0.0;
    spade::io::output_vtk("output", "zer", grid, prim);
    spade::io::binary_read("test.bin", prim);
    spade::io::output_vtk("output", "check", grid, prim);
    
    return 0;
}
