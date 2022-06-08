#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(16, 4, 3);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(32, 32, 48);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    
    cvdf::bound_box_t<double, cvdf::cvdf_dim> bounds;
    
    bounds.min(0) =  0.0;
    bounds.max(0) = 12.0;
    bounds.min(1) = -1.0;
    bounds.max(1) =  1.0;
    bounds.min(2) =  0.0;
    bounds.max(2) =  5.0;
    
    
    cvdf::coords::identity<double> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array flow(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    typedef typename decltype(grid)::dtype real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    real_t p_ref = 570.233265072;
    real_t u_ref = 1202.697;
    real_t t_max = 700.0;
    real_t t_wall = 100.0;
    real_t delta  = 1.0;
    
    real_t alpha = sqrt(1.0 - 100.0/700.0);
    real_t beta = 2.0*alpha*((alpha*alpha-1.0)*atanh(alpha) + alpha)/((alpha*alpha*alpha)*(log(abs(1.0+alpha)) - log(abs(1.0-alpha))));
    
    auto channel_ini = [=](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        cvdf::fluid_state::prim_t<real_t> output(0.0);
        output.p() = p_ref;
        output.u() = u_ref*(1.0 - xyz[1]*xyz[1]/(delta*delta));
        output.v() = 0;
        output.w() = 0;
        output.T() = t_max - (t_max - t_wall)*xyz[1]*xyz[1]/(delta*delta);        
        return output;
    };

    cvdf::algs::fill_array(flow, channel_ini);
    
    bool output = true;
    if (output)
    {
        {
            cvdf::timing::scoped_tmr_t t("export");
            std::string main_filename = cvdf::output::output_vtk("output", "flow", grid, flow);
            if (group.isroot()) print("Exported", main_filename);
        }
    }
    
    grid.exchange_array(flow);
    
    return 0;
}