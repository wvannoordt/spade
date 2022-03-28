#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(3, 3, 3);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(8, 8, 8);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<double, cvdf::cvdf_dim> bounds;
    bounds.min(0) = -1.0;
    bounds.max(0) =  1.0;
    bounds.min(1) = -1.0;
    bounds.max(1) =  1.0;
    bounds.min(2) = -1.0;
    bounds.max(2) =  1.0;
    
    cvdf::coords::identity<double> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array flow(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    typedef typename decltype(grid)::dtype real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    int rank = group.rank();
    auto rank_data = [=](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        cvdf::fluid_state::prim_t<real_t> output(0.0);
        output.p() = rank;
        output.u() = rank;
        output.v() = rank;
        output.w() = rank;
        output.T() = rank;
        return output;
    };

    cvdf::algs::fill_array(flow, rank_data, cvdf::grid::exclude_exchanges);
    
    bool output = false;
    if (output)
    {
        std::string main_filename = cvdf::output::output_vtk("output", "flow", grid, flow);
        if (group.isroot()) print("Exported", main_filename);
    }
    
    grid.exchange_array(flow);
    
    return 0;
}