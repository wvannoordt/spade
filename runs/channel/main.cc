#include "cvdf.h"



int main(int argc, char** argv)
{
    typedef double real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(4, 4, 4);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(56, 56, 56);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    // bounds.min(0) = -delta;
    // bounds.max(0) =  delta;
    // bounds.min(1) =  0.0;
    // bounds.max(1) =  4.0*cvdf::consts::pi*delta;
    // bounds.min(2) =  0.0;
    // bounds.max(2) =  4.0*cvdf::consts::pi*delta/3.0;
    
    // cvdf::coords::integrated_tanh_1D<real_t> yc(bounds.min(1), bounds.max(1), 0.1, 1.3);
    // cvdf::coords::identity_1D<real_t> xc;
    // cvdf::coords::identity_1D<real_t> zc;
    // cvdf::coords::diagonal_coords coords(xc, yc, zc);
    
    cvdf::coords::identity<real_t> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array prim (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array rhs  (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;
    
    // cvdf::convective::totani_lr tscheme(air);
    
    // cvdf::flux_algs::flux_lr_diff(prim, rhs, tscheme);
    
    return 0;
}