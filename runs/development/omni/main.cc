#include <chrono>
#include <spade.h>
#include <PTL.h>

#include "typedef.h"

#include "proto_conv.h"

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
    const int    nxb              = input["Config"]["nx_blck"];
    const int    nyb              = input["Config"]["ny_blck"];
    const int    nguard           = input["Config"]["nguard"];
    const real_t xmin             = input["Config"]["xmin"];
    const real_t xmax             = input["Config"]["xmax"];
    const real_t ymin             = input["Config"]["ymin"];
    const real_t ymax             = input["Config"]["ymax"];
    const bool   do_output        = input["Config"]["output"];
    const std::string init_file   = input["Config"]["init_file"];
    const real_t u0               = input["Fluid"]["u0"];
    const real_t deltau           = input["Fluid"]["deltau"];
    const real_t gamma            = input["Fluid"]["gamma"];
    const real_t b                = input["Fluid"]["b"];
    const real_t cp               = input["Fluid"]["cp"];
    const real_t theta_d          = input["Fluid"]["theta_d"];
    
    //define the gas model
    spade::fluid_state::ideal_gas_t<real_t> air(gamma, (1.0-(1.0/gamma))*cp);
    
    const real_t xc = 0.5*(xmin+xmax);
    const real_t yc = 0.5*(ymin+ymax);
    
    spade::ctrs::array<int, 2> num_blocks(nxb, nyb);
    spade::ctrs::array<int, 2> cells_in_block(nx, ny);
    spade::ctrs::array<int, 2> exchange_cells(nguard, nguard);
    spade::bound_box_t<real_t, 2> bounds;
    bounds.min(0) =  xmin;
    bounds.max(0) =  xmax;
    bounds.min(1) =  ymin;
    bounds.max(1) =  ymax;
    
    
    //restart directory
    std::filesystem::path out_path("checkpoint");
    if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
    
    
    //cartesian coordinate system
    spade::coords::identity<real_t> coords;
    
    //create grid
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    
    
    //create arrays residing on the grid
    prim_t fill1 = 0.0;
    spade::grid::grid_array prim (grid, fill1);
    
    flux_t fill2 = 0.0;
    spade::grid::grid_array rhs  (grid, fill2);
    
    
    //define the initial condition
    const real_t sintheta = std::sin(theta_d*spade::consts::pi/180.0);
    const real_t costheta = std::cos(theta_d*spade::consts::pi/180.0);
    const real_t u_theta  = u0*costheta;
    const real_t v_theta  = u0*sintheta;
    auto ini = [&](const spade::ctrs::array<real_t, 3> x) -> prim_t
    {
        prim_t output;
        const real_t r      = std::sqrt((x[0]-xc)*(x[0]-xc) + (x[1]-yc)*(x[1]-yc));
        const real_t upmax  = deltau*u0;
        const real_t expfac = std::exp(0.5*(1.0-((r*r)/(b*b))));
        const real_t ur     = (1.0/b)*deltau*u0*r*expfac;
        const real_t rhor   = std::pow(1.0 - 0.5*(air.gamma-1.0)*deltau*u0*deltau*u0*expfac, 1.0/(air.gamma - 1.0));
        const real_t pr     = std::pow(rhor, air.gamma)/air.gamma;
        const real_t theta_loc = std::atan2(x[1], x[0]);
        output.p() = pr;
        output.T() = pr/(rhor*air.R);
        output.u() = u_theta - ur*std::sin(theta_loc);
        output.v() = v_theta + ur*std::cos(theta_loc);
        output.w() = 0.0;
        return output;
    };
    
    //fill the initial condition
    spade::algs::fill_array(prim, ini);
    
    //fill the guards
    grid.exchange_array(prim);
    
    //using the 2nd-order centered KEEP scheme
    proto::totani_lr tscheme(air);
    // spade::convective::totani_lr tscheme(air);
    
    
    // spade::pde_algs::flux_div(prim, rhs, tscheme);
    
    // spade::io::output_vtk("output", "rhs", rhs);
    
    return 0;
}
