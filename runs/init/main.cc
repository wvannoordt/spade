#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> num_blocks(16, 4, 3);
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> cells_in_block(32, 32, 32);
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<double, cvdf::cvdf_dim> bounds;
    bounds.min(0) =  0.0;
    bounds.min(1) = -1.0;
    bounds.min(2) =  0.0;
    bounds.max(0) = 12.0;
    bounds.max(1) =  1.0;
    bounds.max(2) =  5.0;
    
    cvdf::coords::identity<double> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array flow(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    typedef typename decltype(grid)::dtype real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    auto channel_ini = [=](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        cvdf::fluid_state::prim_t<real_t> output(0.0);
        output.p() = xyz[1];
        output.u() = xyz[2];
        output.v() = xyz[0];
        output.w() = xyz[1];
        output.T() = xyz[2];
        return output;
    };

#define HAND_ROLL 0

    {cvdf::timing::scoped_tmr_t t("fill");
#if(HAND_ROLL)
        for (int lb = 0; lb < num_blocks[0]*num_blocks[1]*num_blocks[2]; ++lb)
        {
            for (int k = exchange_cells[2]; k < exchange_cells[2]+cells_in_block[2]; ++k)
            {
                for (int j = exchange_cells[1]; j < exchange_cells[1]+cells_in_block[1]; ++j)
                {
                    for (int i = exchange_cells[0]; i < exchange_cells[0]+cells_in_block[0]; ++i)
                    {
                        auto xyz = grid.cell_coords(i, j, k, lb);
                        // flow(0, i, j, k, lb) = 100.0 + sin(xyz[0]);
                        // flow(1, i, j, k, lb) = cos(xyz[0])*sin(xyz[1]);
                        // flow(2, i, j, k, lb) = cos(xyz[1])*sin(xyz[2]);
                        // flow(3, i, j, k, lb) = 0.0;
                        // flow(4, i, j, k, lb) = 100.0 + sin(xyz[0]);
                        
                        auto p  = channel_ini(xyz);
                        flow(0, i, j, k, lb) = p[0];
                        flow(1, i, j, k, lb) = p[1];
                        flow(2, i, j, k, lb) = p[2];
                        flow(3, i, j, k, lb) = p[3];
                        flow(4, i, j, k, lb) = p[4];
                    }
                }
            }
        }
#else
        cvdf::algs::fill_array(flow, channel_ini);
#endif
    }
    
    // std::ofstream myfile("out.vtk");
    // cvdf::output::output_vtk(myfile, grid, flow);
    
    return 0;
}