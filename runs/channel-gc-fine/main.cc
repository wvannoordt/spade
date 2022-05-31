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
    const real_t t_wall = 0.1;
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
                    for (int nnn = 0; nnn < 2; ++nnn)
                    {
                        v4c i_d(ii[0], j-(nnn+0)*nvec_out[1], ii[1], lb[0]);
                        v4c i_g(ii[0], j+(nnn+1)*nvec_out[1], ii[1], lb[0]);
                        prim_t q_d, q_g;
                        for (auto n: range(0,5)) q_d[n[0]] = prims(n[0], i_d[0], i_d[1], i_d[2], i_d[3]);
                        const auto x_g = grid.get_comp_coords(i_g[0], i_g[1], i_g[2], i_g[3]);
                        const auto x_d = grid.get_comp_coords(i_d[0], i_d[1], i_d[2], i_d[3]);
                        const auto n_g = calc_normal_vector(grid.coord_sys(), x_g, i_g, 1);
                        const auto n_d = calc_normal_vector(grid.coord_sys(), x_d, i_d, 1);
                        q_g.p() =  q_d.p();
                        q_g.u() = -q_d.u();
                        q_g.v() = -q_d.v()*n_d[1]/n_g[1];
                        q_g.w() = -q_d.w();
                        q_g.T() =  t_wall;
                        for (auto n: range(0,5)) prims(n[0], i_g[0], i_g[1], i_g[2], i_g[3]) = q_g[n[0]];
                    }
                }
            }
            ++idc;
        }
    }
}

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
        
    bool init_from_file = false;
    std::string init_filename = "";
    cvdf::cli_args::shortname_args_t args(argc, argv);
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
    
    cvdf::ctrs::array<int, cvdf::cvdf_dim> num_blocks(8, 6, 8);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block(36, 24, 36);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*cvdf::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*cvdf::consts::pi*delta;
    
    const real_t targ_cfl = 0.35;
    const int    nt_max   = 150001;
    const int    nt_skip  = 5000;
    const int    checkpoint_skip  = 5000;
    
    cvdf::coords::identity_1D<real_t> xc;
    cvdf::coords::integrated_tanh_1D<real_t> yc(bounds.min(1), bounds.max(1), 0.1, 1.3);
    cvdf::coords::identity_1D<real_t> zc;
    
    cvdf::coords::diagonal_coords coords(xc, yc, zc);
    
    
    
    std::filesystem::path out_path("checkpoint");
    if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
    
    
    cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    cvdf::ctrs::array<real_t, 3> dx_comp(grid.get_dx(0),grid.get_dx(1),grid.get_dx(2));
    cvdf::ctrs::array<real_t, 3> dx_phys
    (
        xc.map(bounds.min(0)+dx_comp[0]) - xc.map(bounds.min(0)),
        yc.map(bounds.min(1)+dx_comp[1]) - yc.map(bounds.min(1)),
        zc.map(bounds.min(2)+dx_comp[2]) - zc.map(bounds.min(2))
    );
    
    
    cvdf::grid::grid_array prim (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array rhs0 (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array rhs1 (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
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
    
    cvdf::ctrs::array<real_t, 3> dx_plus = dx_phys;
    dx_plus *= (rho*u_tau/mu);
    if (group.isroot()) print("dx_p: ", dx_phys);
    if (group.isroot()) print("dx_c: ", dx_comp);
    if (group.isroot()) print("dx^+: ", dx_plus);
    
    const int nidx = 8;
    std::vector<real_t> r_amp_1(cells_in_block[0]/nidx);
    std::vector<real_t> r_amp_2(cells_in_block[1]/nidx);
    std::vector<real_t> r_amp_3(cells_in_block[2]/nidx);
    std::vector<real_t> r_amp_4(grid.get_partition().get_num_local_blocks());
    
    for (auto& p: r_amp_1) p = 1.0 - 2.0*cvdf::utils::unitary_random();
    for (auto& p: r_amp_2) p = 1.0 - 2.0*cvdf::utils::unitary_random();
    for (auto& p: r_amp_3) p = 1.0 - 2.0*cvdf::utils::unitary_random();
    for (auto& p: r_amp_4) p = 1.0 - 2.0*cvdf::utils::unitary_random();
    
    auto ini = [&](const cvdf::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        const real_t shape = 1.0 - pow(x[1]/delta, 4);
        const real_t turb  = du*u_tau*sin(10.0*cvdf::consts::pi*x[1])*cos(12*x[0])*cos(6*x[2]);
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
    
    cvdf::algs::fill_array(prim, ini);
    
    if (init_from_file)
    {
        cvdf::io::binary_read(init_filename, prim);
    }
    
    cvdf::io::output_vtk("output", "ini", grid, prim);
    cvdf::convective::totani_lr tscheme(air);
    // cvdf::convective::weno_3    wscheme(air);
    cvdf::viscous::visc_lr visc_scheme(visc_law);
    
    struct p2c_t
    {
        const cvdf::fluid_state::perfect_gas_t<real_t>* gas;
        typedef prim_t arg_type;
        p2c_t(const cvdf::fluid_state::perfect_gas_t<real_t>& gas_in) {gas = &gas_in;}
        cons_t operator () (const prim_t& q) const
        {
            cons_t w;
            cvdf::fluid_state::convert_state(q, w, *gas);
            return w;
        }
    } p2c(air);
    
    struct c2p_t
    {
        const cvdf::fluid_state::perfect_gas_t<real_t>* gas;
        typedef cons_t arg_type;
        c2p_t(const cvdf::fluid_state::perfect_gas_t<real_t>& gas_in) {gas = &gas_in;}
        prim_t operator () (const cons_t& w) const
        {
            prim_t q;
            cvdf::fluid_state::convert_state(w, q, *gas);
            return q;
        }
    } c2p(air);
    
    struct get_u_t
    {
        const cvdf::fluid_state::perfect_gas_t<real_t>* gas;
        typedef prim_t arg_type;
        get_u_t(const cvdf::fluid_state::perfect_gas_t<real_t>& gas_in) {gas = &gas_in;}
        real_t operator () (const prim_t& q) const
        {
            return sqrt(gas->gamma*gas->R*q.T()) + sqrt(q.u()*q.u() + q.v()*q.v() + q.w()*q.w());
        }
    } get_u(air);
    
    cvdf::reduce_ops::reduce_max<real_t> max_op;
    const real_t time0 = 0.0;
    
    
    
    const real_t dx = cvdf::utils::min(dx_phys[0], dx_phys[1], dx_phys[2]);
    const real_t umax_ini = cvdf::algs::transform_reduce(prim, get_u, max_op);
    const real_t dt     = targ_cfl*dx/umax_ini;
    
    auto ftrans = [&](auto& q) -> void
    {
        cvdf::algs::transform_inplace(prim, p2c);
    };
    auto itrans = [&](auto& q) -> void
    {
        cvdf::algs::transform_inplace(prim, c2p);
    };
    auto calc_rhs = [&](auto& rhs, auto& q, const auto& t) -> void
    {
        rhs = 0.0;
        grid.exchange_array(q);
        set_channel_noslip(q);
        cvdf::flux_algs::flux_lr_diff(q, rhs, tscheme);
	cvdf::flux_algs::flux_lr_diff(q, rhs, visc_scheme);
        cvdf::algs::transform_inplace(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_ar) -> cvdf::ctrs::array<real_t, 5> 
        {
            cvdf::ctrs::array<real_t, 5> rhs_new = rhs_ar;
            rhs_new[2] += force_term;
            return rhs_new;
        });
    };
    
    cvdf::time_integration::rk2 time_int(prim, rhs0, rhs1, time0, dt, calc_rhs, ftrans, itrans);
    
    std::ofstream myfile("hist.dat");
    for (auto nti: range(0, nt_max))
    {
        int nt = nti[0];
        real_t umax   = cvdf::algs::transform_reduce(prim, get_u, max_op);
        if (group.isroot())
        {
            const real_t cfl = umax*dt/dx;
            print(
                "nt: ", cvdf::utils::pad_str(nt, 15),
                "cfl:", cvdf::utils::pad_str(cfl, 15),
                "u+a:", cvdf::utils::pad_str(umax, 15),
                "dx: ", cvdf::utils::pad_str(dx, 15),
                "dt: ", cvdf::utils::pad_str(dt, 15),
                "ftt:", cvdf::utils::pad_str(20.0*u_tau*time_int.time()/delta, 15)
            );
            myfile << nt << " " << cfl << " " << umax << " " << dx << " " << dt << std::endl;
            myfile.flush();
        }
        if (nt%nt_skip == 0)
        {
            if (group.isroot()) print("Output solution...");
            std::string nstr = cvdf::utils::zfill(nt, 8);
            std::string filename = "prims"+nstr;
            cvdf::io::output_vtk("output", filename, grid, prim);
            if (group.isroot()) print("Done.");
        }
        if (nt%checkpoint_skip == 0)
        {
            if (group.isroot()) print("Output checkpoint...");
            std::string nstr = cvdf::utils::zfill(nt, 8);
            std::string filename = "check"+nstr;
            filename = "checkpoint/"+filename+".bin";
            cvdf::io::binary_write(filename, prim);
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
