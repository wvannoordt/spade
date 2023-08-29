#include <stdio.h>
#include "spade.h"

using real_t = double;
using flux_t = spade::fluid_state::flux_t<real_t>;
using prim_t = spade::fluid_state::prim_t<real_t>;
using cons_t = spade::fluid_state::cons_t<real_t>;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    const int nx      = 64;
    const int ny      = 64;
    const int nz      = 64;
    const int nxb     = 4;
    const int nyb     = 4;
    const int nzb     = 4;
    const int nguard  = 2;
    const real_t xmin = 0.0;
    const real_t xmax = 1.0;
    const real_t ymin = 0.0;
    const real_t ymax = 1.0;
    const real_t zmin = 0.0;
    const real_t zmax = 1.0;
    
    spade::ctrs::array<int, 3> num_blocks(nxb, nyb, nzb);
    spade::ctrs::array<int, 3> cells_in_block(nx, ny, nz);
    spade::ctrs::array<int, 3> exchange_cells(nguard, nguard, nguard);
    
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) =  xmin;
    bounds.max(0) =  xmax;
    bounds.min(1) =  ymin;
    bounds.max(1) =  ymax;
    bounds.min(2) =  zmin;
    bounds.max(2) =  zmax;
    
    spade::coords::identity<real_t> coords;
    
    spade::amr::amr_blocks_t blocks(num_blocks, bounds);
    using refine_t = typename decltype(blocks)::refine_type;
    spade::ctrs::array<bool, 3> periodic = true;
    refine_t ref0 = {true, true, true};
    
    spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, blocks, coords, group);
    auto handle = spade::grid::create_exchange(grid, group, periodic);
    
    prim_t fill1 = 0.0;
    flux_t fill2 = 0.0;

    spade::grid::grid_array prim (grid, fill1, spade::device::best);
    spade::grid::grid_array rhs  (grid, fill2, spade::device::best);

    
    const auto ini = _sp_lambda (const spade::coords::point_t<real_t>& x)
    {
        prim_t output;
        output.p() = x[0];
        output.T() = x[1];
        output.u() = x[0];
        output.v() = x[1];
        output.w() = x[2];
        return output;
    };

    {
        spade::timing::scoped_tmr_t t0("fill");
        spade::algs::fill_array(prim, ini);
    }
    
    spade::io::output_vtk("output", "prim", prim);
    
    return 0;
}
