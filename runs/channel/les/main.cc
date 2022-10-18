#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<real_t, 3> v3d;
typedef spade::ctrs::array<real_t, 3> v5d;
typedef spade::ctrs::array<int,    3> v3i;
typedef spade::ctrs::array<int,    4> v4i;
typedef spade::ctrs::array<spade::grid::cell_t<int>, 4> v4c;
typedef spade::fluid_state::prim_t<real_t> prim_t;
typedef spade::fluid_state::flux_t<real_t> flux_t;
typedef spade::fluid_state::cons_t<real_t> cons_t;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    
    bool init_from_file = false;
    std::string init_filename = "";
    spade::cli_args::shortname_args_t args(argc, argv);
    if (args.has_arg("-init"))
    {
        init_filename = args["-init"];
        if (group.isroot()) print("Initialize from", init_filename);
        init_from_file = true;
        if (!std::filesystem::exists(init_filename))
        {
            if (group.isroot()) print("Cannot find ini file", init_filename);
            abort();
        }
    }
    
    v3i filt(4, 8, 3);
    
    spade::ctrs::array<int, 3> num_blocks(8, 8, 8);
    spade::ctrs::array<int, 3> cells_in_block(48/filt[0], 48/filt[1], 48/filt[2]);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    
    spade::bound_box_t<real_t, 3> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*spade::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*spade::consts::pi*delta;
    
    const real_t targ_cfl = 0.25;
    const int    nt_max   = 150000;
    const int    nt_skip  = 10000;
    const int    checkpoint_skip  = 10000;
    
    spade::coords::identity<real_t> coords;
    
    std::filesystem::path out_path("checkpoint");
    if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
    
    
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    prim_t fill1 = 0.0;
    flux_t fill2 = 0.0;
    spade::grid::grid_array prim (grid, fill1);
    spade::grid::grid_array rhs0 (grid, fill2);
    
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
    
    const int nidx = 2;
    std::vector<real_t> r_amp_1(cells_in_block[0]/nidx);
    std::vector<real_t> r_amp_2(cells_in_block[1]/nidx);
    std::vector<real_t> r_amp_3(cells_in_block[2]/nidx);
    std::vector<real_t> r_amp_4(grid.get_partition().get_num_local_blocks());
    
    for (auto& p: r_amp_1) p = 1.0 - 2.0*spade::utils::unitary_random();
    for (auto& p: r_amp_2) p = 1.0 - 2.0*spade::utils::unitary_random();
    for (auto& p: r_amp_3) p = 1.0 - 2.0*spade::utils::unitary_random();
    for (auto& p: r_amp_4) p = 1.0 - 2.0*spade::utils::unitary_random();
    
    auto ini = [&](const spade::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        const real_t shape = 1.0 - pow(x[1]/delta, 4);
        const real_t turb  = du*u_tau*sin(10.0*spade::consts::pi*x[1])*cos(12*x[0])*cos(6*x[2]);
        prim_t output;
        output.p() = p0;
        output.T() = t0;
        output.u() = (20.0*u_tau + 0.0*turb)*shape;
        output.v() = (0.0        + 0.0*turb)*shape;
        output.w() = (0.0        + 0.0*turb)*shape;
        
        int eff_i = i/nidx;
        int eff_j = j/nidx;
        int eff_k = k/nidx;
        
        const real_t per = du*u_tau*(r_amp_1[eff_i] + r_amp_2[eff_j] + r_amp_3[eff_k] + r_amp_4[lb]);
        output.u() += per*shape;
        output.v() += per*shape;
        output.w() += per*shape;
        
        return output;
    };
    
    spade::algs::fill_array(prim, ini);
    
    if (init_from_file)
    {
        spade::io::binary_read(init_filename, prim);
    }
    
    spade::convective::totani_lr tscheme(air);
    spade::convective::weno_3    wscheme(air);
    spade::viscous::visc_lr visc_scheme(visc_law);
    
    struct p2c_t
    {
        const spade::fluid_state::perfect_gas_t<real_t>* gas;
        typedef prim_t arg_type;
        p2c_t(const spade::fluid_state::perfect_gas_t<real_t>& gas_in) {gas = &gas_in;}
        cons_t operator () (const prim_t& q) const
        {
            cons_t w;
            spade::fluid_state::convert_state(q, w, *gas);
            return w;
        }
    } p2c(air);
    
    struct c2p_t
    {
        const spade::fluid_state::perfect_gas_t<real_t>* gas;
        typedef cons_t arg_type;
        c2p_t(const spade::fluid_state::perfect_gas_t<real_t>& gas_in) {gas = &gas_in;}
        prim_t operator () (const cons_t& w) const
        {
            prim_t q;
            spade::fluid_state::convert_state(w, q, *gas);
            return q;
        }
    } c2p(air);
    
    struct get_u_t
    {
        const spade::fluid_state::perfect_gas_t<real_t>* gas;
        typedef prim_t arg_type;
        get_u_t(const spade::fluid_state::perfect_gas_t<real_t>& gas_in) {gas = &gas_in;}
        real_t operator () (const prim_t& q) const
        {
            return sqrt(gas->gamma*gas->R*q.T()) + sqrt(q.u()*q.u() + q.v()*q.v() + q.w()*q.w());
        }
    } get_u(air);
    
    spade::reduce_ops::reduce_max<real_t> max_op;
    real_t time0 = 0.0;
    
    
    
    const real_t dx = spade::utils::min(grid.get_dx(0), grid.get_dx(1), grid.get_dx(2));
    const real_t umax_ini = spade::algs::transform_reduce(prim, get_u, max_op);
    const real_t dt     = targ_cfl*dx/umax_ini;
    
    auto ftrans = [&](auto& q) -> void
    {
        spade::algs::transform_inplace(prim, p2c);
    };
    auto itrans = [&](auto& q) -> void
    {
        spade::algs::transform_inplace(prim, c2p);
    };
    auto calc_rhs = [&](auto& rhs, auto& q, const auto& t) -> void
    {
        print("NEED TO RE_WRTITE");
        abort();
        rhs = 0.0;
        grid.exchange_array(q);
        // set_channel_noslip(q);
        spade::pde_algs::flux_div(q, rhs, tscheme);
        spade::pde_algs::flux_div(prim, rhs, visc_scheme);
        spade::pde_algs::source_term(prim, rhs, [&]()->v5d{return v5d(0,0,force_term,0,0);});
    };
    
    cons_t transform_state;
    spade::fluid_state::state_transform_t trans(prim, transform_state, air);
    
    spade::time_integration::rk2 time_int(prim, rhs0, time0, dt, calc_rhs, trans);
    
    std::ofstream myfile("hist.dat");
    for (auto nti: range(0, nt_max))
    {
        int nt = nti;
        real_t umax   = spade::algs::transform_reduce(prim, get_u, max_op);
        if (group.isroot())
        {
            const real_t cfl = umax*dt/dx;
            print(
                "nt: ", spade::utils::pad_str(nt, 15),
                "cfl:", spade::utils::pad_str(cfl, 15),
                "u+a:", spade::utils::pad_str(umax, 15),
                "dx: ", spade::utils::pad_str(dx, 15),
                "dt: ", spade::utils::pad_str(dt, 15),
                "ftt:", spade::utils::pad_str(20.0*u_tau*time_int.time()/delta, 15)
            );
            myfile << nt << " " << cfl << " " << umax << " " << dx << " " << dt << std::endl;
            myfile.flush();
        }
        if (nt%nt_skip == 0)
        {
            if (group.isroot()) print("Output solution...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "prims"+nstr;
            spade::io::output_vtk("output", filename, grid, prim);
            if (group.isroot()) print("Done.");
        }
        if (nt%checkpoint_skip == 0)
        {
            if (group.isroot()) print("Output checkpoint...");
            std::string nstr = spade::utils::zfill(nt, 8);
            std::string filename = "check"+nstr;
            filename = "checkpoint/"+filename+".bin";
            spade::io::binary_write(filename, prim);
            if (group.isroot()) print("Done.");
        }
        time_int.advance();
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
