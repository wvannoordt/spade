#include <chrono>
#include <spade.h>
#include <PTL.h>

#include "deprecated/time_integration.h"

#include "typedef.h"

int main(int argc, char** argv)
{
    //initialize MPI
    spade::parallel::mpi_t group(&argc, &argv);
    
    //Get the input file
    std::vector<std::string> args;
    for (auto i: range(0, argc)) args.push_back(std::string(argv[i]));
    if (args.size() < 2)
    {
        if (group.isroot()) print("Please provide an input file name!");
        return 1;
    }
    std::string input_filename = args[1];
    
    //read the input file
    PTL::PropertyTree input;
    input.Read(input_filename);
    
    const real_t targ_cfl         = input["Config"]["cfl"];
    const int    nt_max           = input["Config"]["nt_max"];
    const int    nt_skip          = input["Config"]["nt_skip"];
    const int    checkpoint_skip  = input["Config"]["ck_skip"];
    const int    nx               = input["Config"]["nx_cell"];
    const int    ny               = input["Config"]["ny_cell"];
    const int    nz               = input["Config"]["ny_cell"];
    const int    nxb              = input["Config"]["nx_blck"];
    const int    nyb              = input["Config"]["ny_blck"];
    const int    nzb              = input["Config"]["nz_blck"];
    const int    nguard           = input["Config"]["nguard"];
    const real_t xmin             = input["Config"]["xmin"];
    const real_t xmax             = input["Config"]["xmax"];
    const real_t ymin             = input["Config"]["ymin"];
    const real_t ymax             = input["Config"]["ymax"];
    const real_t zmin             = input["Config"]["zmin"];
    const real_t zmax             = input["Config"]["zmax"];
    const std::string init_file   = input["Config"]["init_file"];
    const bool do_output          = input["Config"]["do_output"];
    const bool output_rhs         = input["Config"]["output_rhs"];    
    
    const real_t u0               = input["Fluid"]["u0"];
    const real_t rho0             = input["Fluid"]["rho0"];
    const real_t p0               = input["Fluid"]["p0"];
    const real_t T0               = input["Fluid"]["T0"];
    const real_t mach             = input["Fluid"]["mach"];
    const real_t reynolds         = input["Fluid"]["reynolds"];
    const real_t prandtl          = input["Fluid"]["prandtl"];
    
    
    const real_t r_gas = p0/(rho0*T0);
    const real_t gamma = (u0/mach)*(u0/mach)/(r_gas*T0);
    
    //define the gas model
    spade::fluid_state::ideal_gas_t<real_t> air(gamma, r_gas);
    
    spade::ctrs::array<int, 3> num_blocks(nxb, nyb, nzb);
    spade::ctrs::array<int, 3> cells_in_block(nx, ny, nz);
    spade::ctrs::array<int, 3> exchange_cells(nguard, nguard, nguard);
    
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) =  xmin;
    bounds.max(0) =  xmax;
    bounds.min(1) =  ymin;
    bounds.max(1) =  ymax;
    bounds.min(2) =  zmin;
    bounds.max(2) =  zmax;
    
    const real_t L = 1.0;
    
    
    //restart directory
    std::filesystem::path out_path("checkpoint");
    if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
    
    
    //cartesian coordinate system
    spade::coords::identity<real_t> coords;
    
    //create grid
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    //create arrays residing on the grid
    prim_t fill1 = 0.0;
    spade::grid::grid_array prim(grid, fill1);
    
    flux_t fill2 = 0.0;
    spade::grid::grid_array rhs  (grid, fill2);
    
    print("gah");
    spade::grid::cell_idx_t ii(0,0,0,1);
    prim(0, ii) = 1.0;
    return 1;
    
    //Bull, J. R., Jameson, A. (2015).
    //Simulation of the Taylor-Green vortex using high-order flux reconstruction schemes.
    //AIAA Journal, 53 (9), 2750–2761. https://doi.org/10.2514/1.J053766
    
    //define the initial condition
    auto ini = [&](const spade::ctrs::array<real_t, 3> x) -> prim_t
    {
        prim_t output;
        output.p() = p0 + rho0*u0*u0*(std::cos(2*x[0]/L) + std::cos(2*x[1]/L))*(std::cos(2*x[2]/L) + 2.0);
        output.T() = T0;
        output.u() =  u0*std::sin(x[0]/L)*std::cos(x[1]/L)*std::cos(x[2]/L);
        output.v() = -u0*std::cos(x[0]/L)*std::sin(x[1]/L)*std::cos(x[2]/L);
        output.w() = 0.0;
        return output;
    };
    //fill the initial condition
    spade::algs::fill_array(prim, ini);
    
    //fill the guards
    grid.exchange_array(prim);
    
    //if a restart file is specified, read the data, fill the array, and fill guards
    if (init_file != "none")
    {
        if (group.isroot()) print("reading...");
        spade::io::binary_read(init_file, prim);
        if (group.isroot()) print("Init done.");
        grid.exchange_array(prim);
    }
    
    //using the 2nd-order centered KEEP scheme
    spade::convective::totani_lr        tscheme(air);
    
    //viscous scheme
    const real_t visc = rho0*u0*L/reynolds;
    spade::viscous_laws::constant_viscosity_t<real_t> visc_law(visc, prandtl);
    spade::viscous::visc_lr  vscheme(visc_law);
    
    if (group.isroot())
    {
        print("p0:      ", p0);
        print("rho0:    ", rho0);
        print("T0:      ", T0);
        print("u0:      ", u0);
        print("visc:    ", visc);
        print("prandtl: ", prandtl);
    }
    
    //define an element-wise kernel that returns the acoustic wavespeed for CFL calculation
    struct get_u_t
    {
        const spade::fluid_state::ideal_gas_t<real_t>* gas;
        typedef prim_t arg_type;
        get_u_t(const spade::fluid_state::ideal_gas_t<real_t>& gas_in) {gas = &gas_in;}
        real_t operator () (const prim_t& q) const
        {
            return sqrt(gas->gamma*gas->R*q.T()) + sqrt(q.u()*q.u() + q.v()*q.v() + q.w()*q.w());
        }
    } get_u(air);
    
    spade::reduce_ops::reduce_max<real_t> max_op;
    spade::reduce_ops::reduce_max<real_t> sum_op;
    
    
    
    //calculate timestep
    real_t time0 = 0.0;
    const real_t dx = spade::utils::min(grid.get_dx(0), grid.get_dx(1), grid.get_dx(2));
    const real_t umax_ini = spade::algs::transform_reduce(prim, get_u, max_op);
    const real_t dt     = targ_cfl*dx/umax_ini;
    
    //define the conservative variable transformation
    cons_t transform_state;
    spade::fluid_state::state_transform_t trans(transform_state, air);
    
    //define the residual calculation
    auto calc_rhs = [&](auto& rhs_in, auto& q, const auto& t) -> void
    {
        rhs_in = 0.0;
        grid.exchange_array(q);
        //compute the convective terms
        spade::pde_algs::flux_div(q, rhs_in, tscheme, vscheme);
    };
    
    //define the time integrator
    spade::time_integration::time_axis_t axis(time0, dt);
    spade::time_integration::rk4_t alg;
    spade::time_integration::integrator_data_t q(prim, rhs, alg);
    spade::time_integration::integrator_t time_int(axis, alg, q, calc_rhs, trans);
    
    
    spade::utils::mtimer_t tmr("advance");
    
    //time loop
    for (auto nt: range(0, nt_max+1))
    {
        //cacluate the maximum wavespeed |u|+a
        const real_t umax   = spade::algs::transform_reduce(time_int.solution(), get_u, max_op);
        
        //print some nice things to the screen
        if (group.isroot())
        {
            const real_t cfl = umax*dt/dx;
            const int pn = 10;
            print(
                "nt: ",  spade::utils::pad_str(nt, pn),
                "cfl:",  spade::utils::pad_str(cfl, pn),
                "umax:", spade::utils::pad_str(umax, pn),
                "dx: ",  spade::utils::pad_str(dx, pn),
                "dt: ",  spade::utils::pad_str(dt, pn)
            );
        }
        
        //output the solution
        if (nt%nt_skip == 0)
        {
            if (group.isroot()) print("Output solution...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "prims"+nstr;
            if (do_output) spade::io::output_vtk("output", filename, time_int.solution());
            std::string filename_rhs = "rhs"+nstr;
            if (output_rhs) spade::io::output_vtk("output", filename_rhs, rhs);
            if (group.isroot()) print("Done.");
        }
        
        //output a restart file if needed
        if (nt%checkpoint_skip == 0)
        {
            if (group.isroot()) print("Output checkpoint...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "check"+nstr;
            filename = "checkpoint/"+filename+".bin";
            if (do_output) spade::io::binary_write(filename, time_int.solution());
            if (group.isroot()) print("Done.");
        }
        
        //advance the solution
        tmr.start("advance");
        time_int.advance();
        tmr.stop("advance");
        
        if (group.isroot()) print(tmr);
        
        //check for solution divergence
        if (std::isnan(umax))
        {
            if (group.isroot())
            {
                print("A tragedy has occurred!");
            }
            group.sync();
            return 155;
        }
    }
    return 0;
}
