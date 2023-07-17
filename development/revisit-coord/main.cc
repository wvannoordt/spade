#include "spade.h"

using real_t = double;
using flux_t = spade::fluid_state::flux_t<real_t>;
using prim_t = spade::fluid_state::prim_t<real_t>;
using cons_t = spade::fluid_state::cons_t<real_t>;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    
    spade::ctrs::array<int, 3> num_blocks(4, 4, 4);
    spade::ctrs::array<int, 3> cells_in_block(16, 16, 16);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) =  0.0;
    bounds.max(0) =  2.0;
    bounds.min(1) =  1.0;
    bounds.max(1) =  3.0;
    bounds.min(2) =  -spade::consts::pi/4;
    bounds.max(2) =   spade::consts::pi/4;
    
    // spade::coords::identity<real_t> coords;
    spade::coords::cyl_coords<real_t> coords;
    
    spade::amr::amr_blocks_t blocks(num_blocks, bounds);
    using refine_t = typename decltype(blocks)::refine_type;
    spade::ctrs::array<bool, 3> periodic = false;
    refine_t ref = true;
    const auto crit = [&](const auto& node)
    {
        const auto ii = node.is_domain_boundary();
        return ii.min(1);
    };
    auto nodes = blocks.select(crit);
    blocks.refine(nodes, periodic, ref, spade::amr::constraints::factor2);
    // auto nodez = blocks.select(crit);
    // blocks.refine(nodez, periodic, ref, spade::amr::constraints::factor2);
    
    spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, blocks, coords, group);
    
    prim_t fill1 = 0.0;
    spade::grid::grid_array prim (grid, fill1);
    
    auto ini = [&](const spade::coords::point_t<real_t>& x)
    {
        prim_t output;
        output.p() = x[0];
        output.T() = x[1];
        output.u() = x[2];
        output.v() = 0.0;
        output.w() = 0.0;
        return output;
    };

    spade::algs::fill_array(prim, ini);
    
    spade::io::output_pvtk("output", "prim", prim);
    
    return 0;
}
