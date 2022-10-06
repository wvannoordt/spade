#include <chrono>
#include "spade.h"

typedef double real_t;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3> cells_in_block(16, 16, 16);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    
    bounds.min(0) = 0.0;
    bounds.max(0) = 1.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 1.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 1.0;
    
    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    real_t fill = 0.0;
    
    spade::grid::face_array prim(grid, fill);
    
    // auto ini = [&](const spade::ctrs::array<real_t, 3> x) -> real_t
    // {
    //     const real_t shape = 1.0 - pow(x[1]/delta, 4);
    //     const real_t turb  = du*u_tau*sin(10.0*spade::consts::pi*x[1])*cos(12*x[0])*cos(6*x[2]);
    //     real_t output;
    //     output.p() = p0;
    //     output.T() = t0;
    //     output.u() = (20.0*u_tau + 0.0*turb)*shape;
    //     output.v() = (0.0        + 0.0*turb)*shape;
    //     output.w() = (0.0        + 0.0*turb)*shape;
    // 
    //     int eff_i = i/nidx;
    //     int eff_j = j/nidx;
    //     int eff_k = k/nidx;
    // 
    //     const real_t per = du*u_tau*(r_amp_1[eff_i] + r_amp_2[eff_j] + r_amp_3[eff_k] + r_amp_4[lb]);
    //     output.u() += per*shape;
    //     output.v() += per*shape;
    //     output.w() += per*shape;
    // 
    //     return output;
    // };
    
    // spade::algs::fill_array(prim, ini);
    return 0;
}
