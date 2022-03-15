#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> num_blocks(2, 2, 2);
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> cells_in_block(16, 16, 16);
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<double, cvdf::cvdf_dim> bounds;
    bounds.min(0) = 0.0;
    bounds.min(1) = 0.0;
    bounds.min(2) = 0.0;
    bounds.max(0) = 1.0;
    bounds.max(1) = 1.0;
    bounds.max(2) = 1.0;
    
    cvdf::coords::identity<double> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array flow(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    std::size_t nblocks = num_blocks[0]*num_blocks[1]*num_blocks[2];
    for (auto i: range(0,5,0,cells_in_block[0],0,cells_in_block[1],0,cells_in_block[2],0,nblocks))
    {
        flow(i[0], i[1], i[2], i[3], i[4]) = i[0]+i[1]*i[2]-i[3]-i[4];
    }
    
    std::ofstream myfile("out.vtk");
    cvdf::output::output_grid(myfile, grid, flow);
    
    return 0;
}