#include "cvdf.h"

typedef double real_t;
typedef cvdf::ctrs::array<real_t, 3> v3d;
typedef cvdf::ctrs::array<int,    3> v3i;
typedef cvdf::ctrs::array<int,    4> v4i;
typedef cvdf::ctrs::array<cvdf::grid::cell_t<int>, 4> v4c;
typedef cvdf::fluid_state::prim_t<real_t> prim_t;
typedef cvdf::fluid_state::cons_t<real_t> cons_t;

void set_channel_noslip(auto& prims)
{
    const auto& grid = prims.get_grid();
    for (auto lb: range(0, grid.get_num_local_blocks()))
    {
        const auto& lb_glob = grid.get_partition().get_global_block(lb[0]);
        int idc = 0;
        for (int dir = 2; dir <= 3; ++dir)
        {
            if (grid.is_domain_boundary(lb_glob, dir))
            {
                const auto lb_idx = cvdf::ctrs::expand_index(lb_glob, grid.get_num_blocks());
                const auto nvec_out = v3i(0,2*idc-1,0);
                const cvdf::grid::cell_t<int> j = idc*(grid.get_num_cells(1)-1);
                auto r1 = range(-grid.get_num_exchange(0), grid.get_num_cells(0) + grid.get_num_exchange(0));
                auto r2 = range(-grid.get_num_exchange(2), grid.get_num_cells(2) + grid.get_num_exchange(2));
                for (auto ii: r1*r2)
                {
                    v4c i_domain(ii[0], j,             ii[1], lb[0]);
                    v4c i_ghost (ii[0], j+nvec_out[1], ii[1], lb[0]);
                }
            }
            ++idc;
        }
    }
}

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::ctrs::array<int, cvdf::cvdf_dim> num_blocks(4, 4, 4);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block(24, 24, 24);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*cvdf::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  4.0*cvdf::consts::pi*delta/3.0;
    
    cvdf::coords::identity<real_t> coords;
    
    
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::grid::grid_array prim (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array rhs  (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;
    
    const double p0 = 101325.0;
    const double t0 = 325.0;
    const double u0 = 69.54;
    auto ini = [&](const cvdf::ctrs::array<real_t, 3> x) -> prim_t
    {
        prim_t output;
        output.p() = p0;
        output.T() = t0;
        output.u() = 1.5*(1.0 - x[1]*x[1]/(delta*delta))*u0;
        output.v() = 0.0;
        output.w() = 0.0;
        return output;
    };
    cvdf::algs::fill_array(prim, ini);
    
    cvdf::output::output_vtk("output", "ini", grid, prim);
    
    cvdf::convective::totani_lr tscheme(air);
    cvdf::viscous::visc_lr visc_scheme(visc_law);
    
    const int    nt_max = 1000;
    const real_t dt     = 1e-4;
    for (auto nti: range(0, nt_max))
    {
        int nt = nti[0];
        grid.exchange_array(prim);
        set_channel_noslip(prim);
        cvdf::flux_algs::flux_lr_diff(prim, rhs, tscheme);
        cvdf::flux_algs::flux_lr_diff(prim, rhs, visc_scheme);
        cvdf::algs::transform_inplace(prim, p2c);
        prim += dt*rhs;
        cvdf::algs::transform_inplace(prim, c2p);
        if (group.isroot())
        {
            print(nt);
        }
    }
    
    
    return 0;
}