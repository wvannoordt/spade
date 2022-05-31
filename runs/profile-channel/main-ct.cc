#include "cvdf.h"
#include "local_types.h"
#include "dns_filter.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    
    cvdf::ctrs::array<int, cvdf::cvdf_dim> num_blocks(8, 8, 8);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block(48, 48, 48);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells_filt(8, 8, 8);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*cvdf::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*cvdf::consts::pi*delta;
    
    cvdf::coords::identity<real_t> coords;
    
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    cvdf::grid::cartesian_grid_t grid_filt(num_blocks, cells_in_block, exchange_cells_filt, bounds, coords, group);
    
    cvdf::grid::grid_array prim(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array prim_filt_r(grid_filt, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array prim_filt_f(grid_filt, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;
    
    const real_t p0         = 30.0;
    const real_t t0         = 0.1;
    const real_t u0         = 69.54;
    const real_t mu         = visc_law.get_visc();
    const real_t rho        = p0/(air.R*t0);
    const real_t u_tau      = re_tau*mu/(rho*delta);
    const real_t force_term = rho*u_tau*u_tau/delta;
    const real_t du         = 3.0;
    
    std::vector<std::string> names;
    for (int i = 1; i < argc; i++) names.push_back(std::string(argv[i]));
    std::vector<real_t> y, u;
    int ct = 0;
    for (auto& p: names)
    {
        if (group.isroot()) print(p);
        if (!std::filesystem::exists(p))
        {
            if (group.isroot()) print("The following file does not exsist:", p);
            group.sync();
            return 15;
        }
        cvdf::io::binary_read(p, prim);
        postprocessing::copy_field(prim, prim_filt_r);
        grid_filt.exchange_array(prim_filt_r);
        postprocessing::dns_filter(prim_filt_r, prim_filt_f);
        std::vector<real_t> y_loc, u_loc;
        postprocessing::extract_vel_profile(prim_filt_f, y_loc, u_loc);
        if (y.size()==0)
        {
            y.resize(y_loc.size(), 0.0);
            u.resize(y_loc.size(), 0.0);
        }
        for (int iii = 0; iii < y.size(); ++iii)
        {
            y[iii] += y_loc[iii];
            u[iii] += u_loc[iii];
        }
        ++ct;
    }
    for (auto& v: y) v /= ct;
    for (auto& v: u) v /= ct;
    
    if (group.isroot())
    {
        std::ofstream myfile("u.dat");
        for (int i = 0; i < y.size()/2; ++i)
        {
            myfile << y[i] << " " << u[i] << "\n";
        }
    }
    return 0;
}
