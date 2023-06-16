#include "spade.h"

int main(int argc, char** argv)
{
    using real_t = double;
    
    spade::parallel::mpi_t group(&argc, &argv);
    spade::bound_box_t<real_t, 2> bounds;
    
    bounds.min(0) =  0.0;
    bounds.max(0) =  1.0;
    bounds.min(1) =  0.0;
    bounds.max(1) =  1.0;
    
    spade::coords::identity<real_t> coords;
    spade::ctrs::array<int, 2> num_blocks      (2,  2);
    spade::ctrs::array<int, 2> cells_in_block  (16, 16);
    spade::ctrs::array<int, 2> exchange_cells  (2,  2);
    
    spade::grid::cartesian_blocks_t blocks     (num_blocks, bounds);
    spade::amr::amr_blocks_t        amr_blocks (num_blocks, bounds);
    spade::ctrs::array<bool, 2> per = {true, true};
    amr_blocks.refine(0, per, {true, true});
    amr_blocks.refine(0, per, {true, true});
    
    // const auto get_block_at = [&](const spade::ctrs::array<real_t, 2>& x)
    // {
    //     for (const auto& node: amr_blocks.enumerated_nodes)
    //     {
    //         const auto ub = node->ubox();
    //         if ((ub.min(0) <= x[0] && x[0] < ub.max(0)) && (ub.min(1) <= x[1] && x[1] < ub.max(1)))
    //         {
    //             return node;
    //         }
    //     }
    //     return amr_blocks.enumerated_nodes[0];
    // };
    
    // auto& node = *get_block_at({0.75, 0.5});
    
    spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, amr_blocks, coords, group);
    
    using prim_t = spade::fluid_state::prim_t<real_t>;
    prim_t ft = 0.0;
    spade::grid::cell_array p0(grid, ft);
    spade::grid::cell_array p1(grid, ft);
    
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
    
    spade::algs::fill_array(p0, fill, spade::grid::include_exchanges);
    spade::algs::fill_array(p1, fill, spade::grid::exclude_exchanges);
    
    spade::io::output_vtk("output", "p0", p0);
    spade::io::output_vtk("output", "p1", p1);
    
    return 0;
}