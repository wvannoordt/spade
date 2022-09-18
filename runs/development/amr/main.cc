#include <chrono>
#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<int, 3> v3i;
typedef spade::ctrs::array<int, 2> v2i;
typedef spade::ctrs::array<bool, 3> v3b;
typedef spade::ctrs::array<bool, 2> v2b;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    // spade::bound_box_t<real_t, 2> bounds;
    spade::bound_box_t<real_t, 3> bounds;
    
    const real_t delta = 1.0;
    bounds.min(0) = -delta;
    bounds.max(0) =  delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) = -delta;
    bounds.max(2) =  delta;
    
    // v2i num_blocks(2, 2);
    // v2i cells_in_block(16, 16);
    // v2i exchange_cells(1, 1);
    
    v3i num_blocks(2, 2, 2);
    v3i cells_in_block(16, 16, 16);
    v3i exchange_cells(1, 1, 1);
    
    // spade::block_config::cartesian_blocks_t blocks(num_blocks, bounds);
    spade::amr::amr_blocks_t blocks(num_blocks, bounds);
    
    // v3b refine(true, true, true);
    spade::amr::refine(blocks, *(blocks.all_nodes[0]));
    
    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    spade::io::output_vtk(".", "blocks", blocks, group);
    
    // spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);

    return 0;
}
