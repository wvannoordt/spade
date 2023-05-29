#include "spade.h"

using real_t = double;
using flux_t = spade::fluid_state::flux_t<real_t>;
using prim_t = spade::fluid_state::prim_t<real_t>;
using cons_t = spade::fluid_state::cons_t<real_t>;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    const std::size_t dim = 3;
    spade::ctrs::array<int, dim> num_blocks     (2, 2, 2); 
    spade::ctrs::array<int, dim> cells_in_block (32, 32, 32);
    spade::ctrs::array<int, dim> exchange_cells (2, 2, 2);
    spade::bound_box_t<real_t, dim> bounds;
    for (auto i: range(0, dim))
    {
        bounds.min(i) = 0;
        bounds.max(i) = 1;
    }

    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);

    prim_t fill1 = 0.0;
    flux_t fill2 = 0.0;
    
    spade::grid::grid_array prim(grid, fill1);
    spade::grid::grid_array rhs0(grid, fill2);
    spade::grid::grid_array rhs1(grid, fill2);
    
    spade::fluid_state::ideal_gas_t<real_t> air;
    air.R     = 287.15;
    air.gamma = 1.4;

    using point_type = decltype(grid)::coord_point_type;
    auto ini = [&](const point_type& x) -> prim_t
    {
        prim_t output;
        output.p() = 100.0;
        output.T() = 100.0;
        output.u() = 50 + x[0];
        output.v() = 50 + x[1];
        output.w() = 0.0;
        if constexpr (dim==3)
        {
            output.w() = 50 + x[2];
        }
        
        return output;
    };

    spade::algs::fill_array(prim, ini);
    spade::convective::totani_lr tscheme(air);
    
    spade::bound_box_t<bool, grid.dim()> boundary = false;
    boundary.min(0) = false;
    boundary.max(0) = true;
    boundary.min(1) = false;
    boundary.max(1) = false;

    //Separate boundary and interior calculations
    spade::pde_algs::flux_div     (prim, rhs0, !boundary, tscheme);
    spade::pde_algs::boundary_flux(prim, rhs0,  boundary, tscheme);

    //Single divergence calculation
    spade::pde_algs::flux_div(prim, rhs1, tscheme);

    spade::io::output_vtk("output", "prim", prim);
    spade::io::output_vtk("output", "rhs0", rhs0);
    spade::io::output_vtk("output", "rhs1", rhs1);
    
    return 0;
}
