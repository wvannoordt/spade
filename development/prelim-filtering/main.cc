#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<real_t, 3> v3d;
typedef spade::ctrs::array<int,    3> v3i;
typedef spade::ctrs::array<int,    4> v4i;
typedef spade::ctrs::array<spade::grid::cell_t<int>, 4> v4c;
typedef spade::fluid_state::prim_t<real_t> prim_t;
typedef spade::fluid_state::cons_t<real_t> cons_t;

void test_filter(const auto& src, auto& dest)
{
    const auto& grid = src.get_grid();
    auto rg = grid.get_range(spade::grid::cell_centered);
    const int filtnum = 5;
    spade::ctrs::array<real_t, filtnum> coeff(1.0/filtnum);
    coeff[(filtnum-1)/2] /= 3.0;
    real_t sum = 0.0;
    for (int j = 0; j < filtnum; ++j) sum += coeff[j];
    coeff /= sum;
    coeff /= 3.0;
    for (auto i: range(0,5)*rg)
    {
        dest(i[0], i[1], i[2], i[3], i[4]) = 0.0;
        for (int idir = 0; idir < 3; ++idir)
        {
            spade::ctrs::array<int,4> ijk(i[1], i[2], i[3], i[4]);
            ijk[idir] -= (filtnum-1)/2;
            for (int n = 0; n < filtnum; ++n)
            {
                dest(i[0], i[1], i[2], i[3], i[4]) += src(i[0], ijk[0], ijk[1], ijk[2], ijk[3])*coeff[n];
                ijk[idir]++;
            }
        }
    }
}

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    
    spade::ctrs::array<int, 3> num_blocks(4, 4, 4);
    spade::ctrs::array<int, 3> cells_in_block(48, 48, 48);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    
    const real_t delta = 1.0;
    const real_t re_tau = 180.0;
    
    bounds.min(0) =  0.0;
    bounds.max(0) =  2.0*spade::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  spade::consts::pi*delta;
    
    spade::coords::identity<real_t> coords;   
    
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    spade::grid::grid_array prim_r(grid, 0.0, spade::dims::static_dims<5>(), spade::dims::singleton_dim());
    spade::grid::grid_array prim_f(grid, 0.0, spade::dims::static_dims<5>(), spade::dims::singleton_dim());
    
    spade::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    spade::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;
    
    const real_t p0 = 30.0;
    const real_t t0 = 0.1;
    const real_t u0 = 69.54;
    const real_t mu = visc_law.get_visc();
    const real_t rho = p0/(air.R*t0);
    const real_t u_tau = re_tau*mu/(rho*delta);
    const real_t force_term = rho*u_tau*u_tau/delta;
    const real_t du = 3.0;
    
    std::string prim_filename = "test.bin";
    spade::output::binary_read(prim_filename, prim_r);
    test_filter(prim_r, prim_f);
    spade::output::output_vtk("output", "prim_f", grid, prim_f);
    spade::output::output_vtk("output", "prim_r", grid, prim_r);
    return 0;
}
