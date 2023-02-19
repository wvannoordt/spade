#include <spade.h>

#include "k1.h"

int main(int argc, char** argv)
{
    using real_t = double;
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3> cells_in_block(12, 12, 12);
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
    spade::fluid_state::prim_t<real_t> fill = 0.0;

    spade::grid::grid_array prim (grid, fill);

    local::kernel0_t k0;
    local::kernel1_t k1;

    using stencil_t = local::kernel0_t::stencil_type;
    using array_t   = decltype(prim);
    using data_t    = spade::omni::stencil_data_t<stencil_t, array_t>;

    data_t p;

    return 0;
}
