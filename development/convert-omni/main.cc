#include <iostream>
#include <string>
template <typename thing_t> void print_type(const thing_t& t)
{
    //g++ only
    std::string pf(__PRETTY_FUNCTION__);
    std::size_t start = std::string("void print_type(const thing_t&) [with thing_t = ").length();
    std::size_t end = pf.length()-1;
    std::cout << pf.substr(start, end-start) << std::endl;
}

#include <spade.h>

int main(int argc, char** argv)
{
    using real_t = double;
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3> cells_in_block(10, 10, 10);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) = 0.0;
    bounds.max(0) = 2.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 2.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 2.0;

    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    spade::fluid_state::prim_t<real_t> fill = 0.0;
    spade::grid::grid_array prim (grid, fill);
    spade::grid::grid_array prim2(grid, fill);
    
    using point_type = typename decltype(grid)::coord_point_type;
    using cell_t     = spade::grid::cell_idx_t;

    auto lam = [](const point_type& x, const cell_t& i)
    {
        spade::fluid_state::prim_t<real_t> output;
        output.p() = 1 *x[0] + 2 *x[1] + 3 *x[2] + 10*i.lb();
        output.T() = 4 *x[0] + 5 *x[1] + 6 *x[2];
        output.u() = 7 *x[0] + 8 *x[1] + 9 *x[2];
        output.v() = 10*x[0] + 11*x[1] + 12*x[2];
        output.w() = 13*x[0] + 14*x[1] + 15*x[2];
        return output;
    };

    spade::algs::fill_array(prim, lam);
    spade::io::output_vtk("output", "prim", prim);
    // spade::io::binary_write("test.bin", prim);
    // spade::io::binary_read ("test.bin", prim2);
    // prim2 -= prim;
    // spade::io::output_vtk("output", "diff", prim2);

    return 0;
}
