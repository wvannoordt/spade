#include "symd.h"
#include "spade.h"

int main(int argc, char** argv)
{
    using real_t = double;
    
    spade::parallel::mpi_t group(&argc, &argv);
    spade::bound_box_t<real_t, 2> bounds;
    
    bounds.min(0) = 0.0;
    bounds.max(0) = 1.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 1.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 1.0;
    
    spade::coords::identity<real_t> coords;
    spade::ctrs::array<int, 2> num_blocks      (2,  2);
    spade::ctrs::array<int, 2> cells_in_block  (32, 32);
    spade::ctrs::array<int, 2> exchange_cells  (2,  2);
    
    spade::grid::cartesian_blocks_t blocks     (num_blocks, bounds);
    spade::amr::amr_blocks_t        amr_blocks (num_blocks, bounds);
    amr_blocks.refine(3, {false, false}, {true, true});
    // print(amr_blocks.total_num_blocks());
    
    amr_blocks.all_nodes[7]->debug();
    // print(amr_blocks.total_num_blocks());
    // print(amr_blocks.all_nodes[0]->is_domain_boundary());
    
    spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, amr_blocks, coords, group);
    
    using prim_t = spade::fluid_state::prim_t<real_t>;
    prim_t ft = 0.0;
    spade::grid::cell_array prim(grid, ft);
    
    using point_t = decltype(grid)::coord_point_type;
    auto fill = [](const point_t& x)
    {
        prim_t output;
        output.p() = 10.0;
        output.T() = 10.0;
        output.u() = 10.0*std::sin(2.0*spade::consts::pi*x[0] + 2.0*spade::consts::pi*x[1]);
        output.v() = 0.0;
        output.w() = 0.0;
        return output;
    };
    
    spade::algs::fill_array(prim, fill);
    
    spade::io::output_vtk("output", "prim", prim);
    
    return 0;
}