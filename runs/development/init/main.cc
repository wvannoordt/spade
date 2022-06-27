#include "spade.h"

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    
    spade::ctrs::array<std::size_t, 3> num_blocks(16, 4, 3);
    spade::ctrs::array<std::size_t, 3> cells_in_block(32, 32, 48);
    spade::ctrs::array<std::size_t, 3> exchange_cells(2, 2, 2);
    
    spade::bound_box_t<double, 3> bounds;
    
    bounds.min(0) =  0.0;
    bounds.max(0) = 12.0;
    bounds.min(1) = -1.0;
    bounds.max(1) =  1.0;
    bounds.min(2) =  0.0;
    bounds.max(2) =  5.0;
    
    
    spade::coords::identity<double> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    spade::grid::grid_array flow(grid, 0.0, spade::dims::static_dims<5>(), spade::dims::singleton_dim());
    
    typedef typename decltype(grid)::dtype real_t;
    typedef spade::ctrs::array<real_t, 3> v3d;
    
    real_t p_ref = 570.233265072;
    real_t u_ref = 1202.697;
    real_t t_max = 700.0;
    real_t t_wall = 100.0;
    real_t delta  = 1.0;
    
    real_t alpha = sqrt(1.0 - 100.0/700.0);
    real_t beta = 2.0*alpha*((alpha*alpha-1.0)*atanh(alpha) + alpha)/((alpha*alpha*alpha)*(log(abs(1.0+alpha)) - log(abs(1.0-alpha))));
    
    auto channel_ini = [=](const v3d& xyz) -> spade::fluid_state::prim_t<real_t>
    {
        spade::fluid_state::prim_t<real_t> output(0.0);
        output.p() = p_ref;
        output.u() = u_ref*(1.0 - xyz[1]*xyz[1]/(delta*delta));
        output.v() = 0;
        output.w() = 0;
        output.T() = t_max - (t_max - t_wall)*xyz[1]*xyz[1]/(delta*delta);        
        return output;
    };

    spade::algs::fill_array(flow, channel_ini);
    
    bool output = true;
    if (output)
    {
        {
            spade::timing::scoped_tmr_t t("export");
            std::string main_filename = spade::output::output_vtk("output", "flow", grid, flow);
            if (group.isroot()) print("Exported", main_filename);
        }
    }
    
    grid.exchange_array(flow);
    
    return 0;
}