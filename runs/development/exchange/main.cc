#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(3, 3, 2);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(16, 16, 16);
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
    
    const double pi = 3.14159265359;
    
    int rank = group.rank();
    auto rank_data = [=](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        cvdf::fluid_state::prim_t<real_t> output(0.0);
        double val = sin(2.0*pi*xyz[0])+2.0*cos(8.0*pi*xyz[1])*sin(4.0*pi*xyz[2]);
        output[0] = val;
        output[1] = val;
        output[2] = val;
        output[3] = val;
        output[4] = val;
        return output;
    };

    cvdf::algs::fill_array(flow, rank_data, cvdf::grid::exclude_exchanges);
    grid.exchange_array(flow);
    
    double err = 0.0;
    std::size_t num_ct = 0;
    for (auto i: grid.get_range(flow.centering_type(), cvdf::grid::include_exchanges))
    {
        auto xyz = grid.cell_coords(i[0], i[1], i[2], i[3]);
        double val = sin(2.0*pi*xyz[0])+2.0*cos(8.0*pi*xyz[1])*sin(4.0*pi*xyz[2]);
        double err_loc = flow(0, i[0], i[1], i[2], i[3]) - val;
        err_loc *= err_loc;
        err = err + err_loc;
    }
    err = group.sum(err);
    
    if (group.isroot()) print("Accumulated error:", err);
    
    bool output = false;
    if (output)
    {
        std::string main_filename = cvdf::output::output_vtk("output", "flow", grid, flow);
        if (group.isroot()) print("Exported", main_filename);
    }
    
    
    return 0;
}