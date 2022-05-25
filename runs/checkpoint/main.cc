#include "cvdf.h"

typedef double real_t;
typedef cvdf::ctrs::array<real_t, 3> v3d;
typedef cvdf::ctrs::array<int,    3> v3i;
typedef cvdf::ctrs::array<int,    4> v4i;
typedef cvdf::ctrs::array<cvdf::grid::cell_t<int>, 4> v4c;
typedef cvdf::fluid_state::prim_t<real_t> prim_t;
typedef cvdf::fluid_state::cons_t<real_t> cons_t;

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<int, cvdf::cvdf_dim> num_blocks(2, 2, 2);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block(16, 16, 16);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    bounds.min(0) = 0.0;
    bounds.max(0) = 1.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 1.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 1.0;
    
    cvdf::coords::identity<real_t> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array prim (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    auto ini = [&](const cvdf::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        prim_t output;
        output.p() = 0.5*x[0];
        output.T() = 0.5*x[1];
        output.u() = 0.5*sin(x[0]);
        output.v() = 0.9*cos(x[1]);
        output.w() = 1.2*cos(x[2]);
        return output;
    };
    
    cvdf::algs::fill_array(prim, ini);
    cvdf::output::output_vtk("output", "ini", grid, prim);
    cvdf::output::binary_write("test.bin", prim);
    prim = 0.0;
    cvdf::output::output_vtk("output", "zer", grid, prim);
    cvdf::output::binary_read("test.bin", prim);
    cvdf::output::output_vtk("output", "check", grid, prim);
    
    return 0;
}