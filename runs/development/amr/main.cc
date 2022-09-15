#include <chrono>
#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<int, 3> v3i;
typedef spade::ctrs::array<int, 2> v2i;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    // spade::bound_box_t<real_t, 2> bounds;
    spade::bound_box_t<real_t, 3> bounds;
    
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*spade::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*spade::consts::pi*delta;
    
    // v2i num_blocks(5, 10);
    // v2i cells_in_block(16, 16);
    // v2i exchange_cells(1, 1);
    
    v3i num_blocks(4, 3, 5);
    v3i cells_in_block(16, 16, 16);
    v3i exchange_cells(1, 1, 1);
    
    // spade::block_config::cartesian_blocks_t blocks(num_blocks, bounds);
    spade::block_config::amr_blocks_t blocks(num_blocks, bounds);
    
    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    spade::io::output_vtk(".", "blocks", blocks, group);
    // spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);

    return 0;
}
