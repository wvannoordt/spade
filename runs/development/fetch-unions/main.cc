#include <chrono>
#include "spade.h"

using real_t = double;
using prim_t = spade::fluid_state::prim_t<real_t>;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);    
    spade::ctrs::array<int, 2> num_blocks(2, 2);
    spade::ctrs::array<int, 2> cells_in_block(16, 16);
    spade::ctrs::array<int, 2> exchange_cells(2, 2);
    spade::bound_box_t<real_t, 2> bounds;
    bounds.min(0) = 0.0;
    bounds.max(0) = 1.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 1.0;
    
    spade::coords::identity<real_t> coords;  
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    
    prim_t fill = 0.0;    
    spade::grid::grid_array prim(grid, fill);
    spade::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;

    spade::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;

    auto ini = [&](const spade::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        real_t yy = 2.0*x[1]-1.0;
        prim_t output;
        output.p() = 100.0;
        output.T() = 300.0;
        output.u() = 70*(1.0-yy*yy);
        output.v() = 0.0;
        output.w() = 0.0;
        return output;
    };
    
    spade::algs::fill_array(prim, ini);
    spade::convective::totani_lr tscheme(air);
    spade::convective::weno_3    wscheme(air);
    spade::viscous::visc_lr  visc_scheme(visc_law);
    
    spade::fetch::face_fetch_t
    <
        spade::fetch::flux_line
        <
            4,
            spade::fetch::cell_info
            <
                spade::fetch::cell_state<spade::fluid_state::prim_t<real_t>>,
                spade::fetch::cell_normal<spade::ctrs::array<real_t, 3>>
            >
        >,
        spade::fetch::face_info<>
    > info;
    
    spade::grid::cell_idx_t icell(4, 4, 0, 0);
    spade::grid::face_idx_t iface = spade::grid::cell_to_face(icell, 1, 1);
    
    spade::fetch::get_flux_data(grid, prim, iface, info);
    print(info);
    
    return 0;
}