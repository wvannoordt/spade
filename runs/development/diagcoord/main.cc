#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(6, 3, 2);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(32, 32, 32);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<double, cvdf::cvdf_dim> bounds;
    bounds.min(0) =  0.0;
    bounds.max(0) =  8.0;
    bounds.min(1) = -1.0;
    bounds.max(1) =  1.0;
    bounds.min(2) =  0.0;
    bounds.max(2) =  4.0;
    
    const double eta1 = bounds.max(1);
    const double eta0 = bounds.min(1);
    const double deta = eta1 - eta0;
    const double k = 1.3;
    const double strch = 0.1;
    const double alpha0 = k/deta;
    const double alpha1 = -k/deta;
    const double beta0  =  alpha0*(strch*deta-eta0);
    const double beta1  =  alpha0*(eta1+strch*deta);
    
    auto abs = [](const double& d) -> double {return d<0?-d:d;};
    auto func = [=](const double& eta) -> double
    {
        return  log(abs(cosh(alpha0*eta+beta0)))/alpha0  +  log(abs(cosh(alpha1*eta+beta1)))/alpha1 - eta;
    };
    
    const double f0 = func(eta0);
    const double norm = func(eta1)-f0;
    
    auto mp = [=](const double& x) -> double {return eta0 + deta*(func(x) - f0)/norm;};
    cvdf::coords::analytical_1D yc(mp);
    cvdf::coords::identity_1D<double> xc;
    cvdf::coords::identity_1D<double> zc;
    cvdf::coords::diagonal_coords coords(xc, yc, zc);
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array cons(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    typedef typename decltype(grid)::dtype real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    auto channel_ini = [=](const v3d& xyz) -> cvdf::fluid_state::cons_t<real_t>
    {
        cvdf::fluid_state::cons_t<real_t> output(0.0);
        return output;
    };

    cvdf::algs::fill_array(cons, channel_ini, cvdf::grid::include_exchanges);
    
    bool output = false;
    if (output)
    {
        std::string main_filename = cvdf::output::output_vtk("output", "ini", grid, cons);
        if (group.isroot()) print("Exported", main_filename);
    }
    
    // cvdf::time_integration::time_integrator(flow)
    // cvdf::navier_stokes::convective_operator<>
    
    return 0;
}