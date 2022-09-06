#include <chrono>
#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<int,    3> v3i;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
        
    v3i num_blocks(8, 8, 8);
    v3i cells_in_block(16, 16, 16);
    v3i exchange_cells(1, 1, 1);
    spade::bound_box_t<real_t, 3> bounds;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*spade::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*spade::consts::pi*delta;

    spade::coords::identity<real_t> coords;
    spade::block_config::cartesian_blocks_t blocks(num_blocks, bounds);
    // spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, blocks, bounds, coords, group);

    return 0;
}
