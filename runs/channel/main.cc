#include "cvdf.h"

int main(int argc, char** argv)
{
    typedef double real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    // cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(6, 3, 2);
    // cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(32, 32, 32);
    // cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    // cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    // bounds.min(0) =  0.0;
    // bounds.max(0) =  8.0;
    // bounds.min(1) = -1.0;
    // bounds.max(1) =  1.0;
    // bounds.min(2) =  0.0;
    // bounds.max(2) =  4.0;
    
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(4, 4, 4);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(56, 56, 56);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    
    bounds.min(0) = -1.0;
    bounds.max(0) =  1.0;
    bounds.min(1) = -1.0;
    bounds.max(1) =  1.0;
    bounds.min(2) = -1.0;
    bounds.max(2) =  1.0;
    
    // cvdf::coords::integrated_tanh_1D<real_t> yc(bounds.min(1), bounds.max(1), 0.1, 1.3);
    // cvdf::coords::identity_1D<real_t> xc;
    // cvdf::coords::identity_1D<real_t> zc;
    // cvdf::coords::diagonal_coords coords(xc, yc, zc);
    
    cvdf::coords::identity<real_t> coords;
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array prim    (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array rhs     (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array rhs_test(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;
    
    real_t p_ref = 570.233265072;
    real_t u_ref = 1202.697;
    real_t t_max = 700.0;
    real_t t_wall = 100.0;
    real_t delta  = 1.0;
    
    real_t alpha = sqrt(1.0 - 100.0/700.0);
    real_t beta = 2.0*alpha*((alpha*alpha-1.0)*atanh(alpha) + alpha)/((alpha*alpha*alpha)*(log(abs(1.0+alpha)) - log(abs(1.0-alpha))));
    
    auto channel_ini = [=](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        cvdf::fluid_state::prim_t<real_t> temp;
        temp.p() = p_ref;
        temp.u() = u_ref*(1.0 - xyz[1]*xyz[1]/(delta*delta));
        temp.v() = 0;
        temp.w() = 0;
        temp.T() = t_max - (t_max - t_wall)*xyz[1]*xyz[1]/(delta*delta);
        return temp;
    };
    
    cvdf::navier_stokes_mms::cns_pergectgas_mms_t mms(air, visc_law);
    auto mms_test_func = [&](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        return mms.test_fcn(xyz);
    };
    
    auto mms_conv_func = [&](const v3d& xyz) -> cvdf::ctrs::array<real_t,5>
    {
        return mms.conv_rhs(xyz);
    };

    cvdf::algs::fill_array(prim,     mms_test_func, cvdf::grid::include_exchanges);
    cvdf::algs::fill_array(rhs_test, mms_conv_func, cvdf::grid::include_exchanges);
    
    
    // std::string main_filename = cvdf::output::output_vtk("output", "ini", grid, prim);
    // if (group.isroot()) print("Exported", main_filename);    
    
    cvdf::convective::totani_lr tscheme(air);
    cvdf::flux_algs::flux_lr_diff(prim, rhs, tscheme);
    cvdf::output::output_vtk("output", "rhs_num", grid, rhs);
    cvdf::output::output_vtk("output", "rhs_ana", grid, rhs_test);
    rhs -= rhs_test;
    cvdf::output::output_vtk("output", "rhs_err", grid, rhs);
    
    return 0;
}