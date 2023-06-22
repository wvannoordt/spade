
#include "spade.h"

int main(int argc, char** argv)
{
    using real_t = double;
    
    spade::parallel::mpi_t group(&argc, &argv);
    
    spade::bound_box_t<real_t, 2> bounds;
    
    bounds.min(0) =   0.0;
    bounds.max(0) =   1.0;
    bounds.min(1) =   0.0;
    bounds.max(1) =   1.0;
    
    spade::coords::identity<real_t> coords;
    spade::ctrs::array<int, 2> num_blocks      (4,  4);
    spade::ctrs::array<int, 2> cells_in_block  (16, 16);
    spade::ctrs::array<int, 2> exchange_cells  (2,  2);
    
    spade::grid::cartesian_blocks_t blocks     (num_blocks, bounds);
    spade::amr::amr_blocks_t        amr_blocks (num_blocks, bounds);
    spade::ctrs::array<bool, 2> per = {false, false};
    
    using ref_t = typename decltype(amr_blocks)::refine_type;
    ref_t ref = {true, true};
    
    const auto constraint = spade::amr::constraints::factor2;
    
    for (int i = 0; i < 3; ++i)
    {
        const auto selection = [&](const auto& node)
        {
            const auto& bbox = amr_blocks.get_block_box(node.tag);
            spade::ctrs::array<real_t, 2> c0(bbox.min(0), bbox.min(1));
            spade::ctrs::array<real_t, 2> c1(bbox.max(0), bbox.min(1));
            spade::ctrs::array<real_t, 2> c2(bbox.min(0), bbox.max(1));
            spade::ctrs::array<real_t, 2> c3(bbox.max(0), bbox.max(1));
            
            const auto in = [](const auto& c) {return spade::ctrs::array_norm(c) < 0.35;};
            
            bool b0 = in(c0);
            bool b1 = in(c1);
            bool b2 = in(c2);
            bool b3 = in(c3);
            
            return !((b0&&b1&&b2&&b3) || (!b0&&!b1&&!b2&&!b3));
        };
        auto nodes = amr_blocks.select(selection);
        amr_blocks.refine(nodes, per, ref, constraint);
    }
    
    spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, amr_blocks, coords, group);
    
    auto handle = spade::grid::create_exchange(grid, group, per);
    
    using prim_t = spade::fluid_state::prim_t<real_t>;
    prim_t ft = 0.0;
    spade::grid::cell_array p0(grid, ft);
    spade::grid::cell_array p1(grid, ft);
    
    using point_t = decltype(grid)::coord_point_type;
    auto fill = [](const point_t& x, const spade::grid::cell_idx_t ii)
    {
        prim_t output;
        output.p() = -10.0;
        output.T() = -10.0;
        output.u() = 10.0*std::sin(2.0*spade::consts::pi*x[0] + 2.0*spade::consts::pi*x[1]);
        output.v() = ii.lb();
        output.w() = -10.0;
        return output;
    };
    
    spade::algs::fill_array(p0, fill, spade::grid::include_exchanges);
    spade::algs::fill_array(p1, fill, spade::grid::exclude_exchanges);
    
    spade::io::output_vtk("output", "p0", p0);
    spade::io::output_vtk("output", "p1", p1);
    
    handle.exchange(p1);
    
    spade::io::output_vtk("output", "exch", p1);
    return 0;
}