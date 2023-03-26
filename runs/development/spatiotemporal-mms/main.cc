#include <chrono>
#include "spade.h"

using real_t = double;

struct sadv_t
{
    using omni_type = spade::omni::stencil_t<
        spade::grid::face_centered,
        spade::omni::elem_t<
            spade::omni::offset_t<-3,0,0>,
            spade::omni::info_list_t<spade::omni::info::value>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t<-1,0,0>,
            spade::omni::info_list_t<spade::omni::info::value>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t< 1,0,0>,
            spade::omni::info_list_t<spade::omni::info::value>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t< 3,0,0>,
            spade::omni::info_list_t<spade::omni::info::value>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t<0,0,0>,
            spade::omni::info_list_t<spade::omni::info::normal>
        >
    >;

    spade::ctrs::array<real_t, 2> velocity;

    real_t operator() (const auto& input) const
    {
        const auto& n  = spade::omni::access<spade::omni::info::normal>(input.face(0_c));
        const auto& p0 = spade::omni::access<spade::omni::info::value >(input.cell(0_c));
        const auto& p1 = spade::omni::access<spade::omni::info::value >(input.cell(1_c));
        const auto& p2 = spade::omni::access<spade::omni::info::value >(input.cell(2_c));
        const auto& p3 = spade::omni::access<spade::omni::info::value >(input.cell(3_c));
        const auto& u  = velocity[0]*n[0] + velocity[1]*n[1];
        const real_t c0 = 0.0;
        const real_t c1 = 0.5;
        const real_t c2 = 0.5;
        const real_t c3 = 0.0;
        return c0*p0 + c1*p1 + c2*p2 + c3*p3;
    }
};

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);

    spade::ctrs::array<int, 2> num_blocks(4, 4);
    spade::ctrs::array<int, 2> cells_in_block(64, 64);
    spade::ctrs::array<int, 2> exchange_cells(2, 2);
    spade::bound_box_t<real_t, 2> bounds;
    bounds.min(0) =  -1.0;
    bounds.max(0) =   1.0;
    bounds.min(1) =  -1.0;
    bounds.max(1) =   1.0;
    
    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    
    real_t fill1 = 0.0;

    spade::grid::grid_array phi (grid, fill1);
    spade::grid::grid_array rhs (grid, fill1);


    auto ini = [&](const spade::coords::point_t<real_t>& x)
    {
        const real_t r2 = x[1]*x[1]+x[0]*x[0];
        return std::exp(-5.0*r2);
    };
    
    spade::algs::fill_array(phi, ini);
    grid.exchange_array(phi);
    
    spade::ctrs::array<real_t, 2> vel(0.5, 0.0);
    const real_t targ_cfl = 0.6;
    const real_t dx       = spade::utils::min(grid.get_dx(0), grid.get_dx(1), grid.get_dx(2));
    const real_t umax_ini = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
    const real_t dt       = targ_cfl*dx/umax_ini;


    sadv_t sadv{vel};
    auto calc_rhs = [&](auto& resid, const auto& sol, const auto& t) -> void
    {
        resid = 0.0;
        spade::pde_algs::flux_div(sol, resid, sadv);
    };
    
    auto boundary_cond = [&](auto& sol, const auto& t) -> void
    {
        grid.exchange_array(sol);
    };

    real_t time0 = 0.0;
    spade::time_integration::time_axis_t axis(time0, dt);
    spade::time_integration::ssprk34_t alg;
    spade::time_integration::integrator_data_t q(phi, rhs, alg);
    spade::time_integration::integrator_t time_int(axis, alg, q, calc_rhs, boundary_cond);
    spade::timing::mtimer_t tmr("advance");
    const int nt_max  = 10000;
    const int nt_skip = nt_max/100;
    for (auto nt: range(0, nt_max+1))
    {
        if (nt%nt_skip == 0)
        {
            if (group.isroot()) print("Output solution...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "prims"+nstr;
            spade::io::output_vtk("output", filename, time_int.solution());
            if (group.isroot()) print("Done.");
        }
    	tmr.start("advance");
        time_int.advance();
        tmr.stop("advance");
        if (group.isroot()) print(tmr);
    }
    
    return 0;
}
