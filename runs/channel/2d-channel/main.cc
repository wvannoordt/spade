#include "spade.h"
#include "get_profile.h"

const int dim = 2;

void set_channel_noslip(auto& prims)
{
    const real_t t_wall = 0.1;
    const auto& grid = prims.get_grid();
    for (auto lb: range(0, grid.get_num_local_blocks()))
    {
        const auto& lb_glob = grid.get_partition().get_global_block(lb);
        int idc = 0;
        for (int dir = 2; dir <= 3; ++dir)
        {
            if (grid.is_domain_boundary(lb_glob, dir))
            {
                const auto lb_idx = spade::ctrs::expand_index(lb_glob, grid.get_num_blocks());
                const auto nvec_out = v3i(0,2*idc-1,0);
                const spade::grid::cell_t<int> j = idc*(grid.get_num_cells(1)-1);
                auto r1 = range(-grid.get_num_exchange(0), grid.get_num_cells(0) + grid.get_num_exchange(0));
                auto r2 = range(-grid.get_num_exchange(2), grid.get_num_cells(2) + grid.get_num_exchange(2));
                for (auto ii: r1*r2)
                {
                    for (int nnn = 0; nnn < 2; ++nnn)
                    {
                        v4c i_d(ii[0], j-(nnn+0)*nvec_out[1], ii[1], lb);
                        v4c i_g(ii[0], j+(nnn+1)*nvec_out[1], ii[1], lb);
                        prim_t q_d, q_g;
                        for (auto n: range(0,5)) q_d[n] = prims(n, i_d[0], i_d[1], i_d[2], i_d[3]);
                        const auto x_g = grid.get_comp_coords(i_g[0], i_g[1], i_g[2], i_g[3]);
                        const auto x_d = grid.get_comp_coords(i_d[0], i_d[1], i_d[2], i_d[3]);
                        const auto n_g = calc_normal_vector(grid.coord_sys(), x_g, i_g, 1);
                        const auto n_d = calc_normal_vector(grid.coord_sys(), x_d, i_d, 1);
                        q_g.p() =  q_d.p();
                        q_g.u() = -q_d.u();
                        q_g.v() = -q_d.v()*n_d[1]/n_g[1];
                        q_g.w() = -q_d.w();
                        q_g.T() =  t_wall;
                        for (auto n: range(0,5)) prims(n, i_g[0], i_g[1], i_g[2], i_g[3]) = q_g[n];
                    }
                }
            }
            ++idc;
        }
    }
}

void garbo_visc(auto& prims, auto& rhs, real_t mu)
{
    auto& grd = prims.get_grid();
    auto rg = grd.get_range(spade::grid::cell_centered);
    real_t dx = grd.get_dx(0);
    real_t dy = grd.get_dx(1);
    real_t dz = grd.get_dx(2);
    for (auto idir: range(0, 3))
    {
        for (auto i: rg)
        {
            real_t lap = 0.0;
            real_t d00 = prims(2+idir, i[0],   i[1],   i[2],   i[3]);
            real_t d0p = prims(2+idir, i[0]+1, i[1],   i[2],   i[3]);
            real_t d0m = prims(2+idir, i[0]-1, i[1],   i[2],   i[3]);
            real_t d1p = prims(2+idir, i[0],   i[1]+1, i[2],   i[3]);
            real_t d1m = prims(2+idir, i[0],   i[1]-1, i[2],   i[3]);
            if constexpr (dim==3)
            {
                real_t d2p = prims(2+idir, i[0],   i[1],   i[2]+1, i[3]);
                real_t d2m = prims(2+idir, i[0],   i[1],   i[2]-1, i[3]);
                lap += (d2p - 2*d00 + d2m)/(dz*dz);
            }
            lap += (d0p - 2*d00 + d0m)/(dx*dx);
            lap += (d1p - 2*d00 + d1m)/(dy*dy);
            lap *= mu;
            rhs(2+idir, i[0],   i[1],   i[2],   i[3]) += lap;
        }
    }
}

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
    
    spade::ctrs::array<int, dim> num_blocks(4,2);
    spade::ctrs::array<int, dim> cells_in_block(48, 48);
    spade::ctrs::array<int, dim> exchange_cells(2, 2);
    //spade::ctrs::array<int, 2> exchange_cells(8, 8, 8);
    spade::bound_box_t<real_t, dim> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*spade::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*spade::consts::pi*delta;
    
    const real_t targ_cfl = 0.45;
    const int    nt_max   = 30000001;
    const int    nt_skip  = 10000;
    const int    checkpoint_skip  = 10000;
    
    
    // spade::coords::identity<real_t> coords;
    
    spade::coords::identity_1D<real_t> xc;
    spade::coords::integrated_tanh_1D<real_t> yc(bounds.min(1), bounds.max(1), 0.1, 1.3);
    spade::coords::identity_1D<real_t> zc;
    spade::coords::diagonal_coords coords(xc, yc, zc);
    
    std::filesystem::path out_path("checkpoint");
    if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group, spade::static_math::int_const_t<dim>());
    
    
    prim_t fill = 0.0;
    spade::grid::grid_array prim (grid, fill);
    spade::grid::grid_array rhs0 (grid, fill);
    spade::grid::grid_array rhs1 (grid, fill);    
    
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
    const real_t mach = 20.0*u_tau/sqrt(air.R*air.gamma*t0);
    
    const int nidx = 8;
    std::vector<real_t> r_amp_1(spade::utils::max(cells_in_block[2]/nidx, 10));
    std::vector<real_t> r_amp_2(spade::utils::max(cells_in_block[2]/nidx, 10));
    std::vector<real_t> r_amp_3(spade::utils::max(cells_in_block[2]/nidx, 10));
    std::vector<real_t> r_amp_4(grid.get_partition().get_num_local_blocks());
    
    for (auto& p: r_amp_1) p = 1.0 - 2.0*spade::utils::unitary_random();
    for (auto& p: r_amp_2) p = 1.0 - 2.0*spade::utils::unitary_random();
    for (auto& p: r_amp_3) p = 1.0 - 2.0*spade::utils::unitary_random();
    for (auto& p: r_amp_4) p = 1.0 - 2.0*spade::utils::unitary_random();
    int i3d = ((dim==3)?1:0);
    auto ini = [&](const spade::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        prim_t output;
        output.p() = p0;
        output.T() = t0;
        output.u() = 0.5*re_tau*u_tau*(1.0-(x[1]*x[1])/(delta*delta));
        output.v() = 0.0;
        output.w() = 0.0;
        return output;
    };
    
    spade::algs::fill_array(prim, ini);
    if (init_from_file)
    {
        spade::io::binary_read(init_filename, prim);
        std::vector<real_t> yprof, uprof;
        get_profile(prim, yprof, uprof);
        if (group.isroot())
        {
            std::ofstream myfile("prof.dat");
            for (int ppp = 0; ppp < yprof.size(); ++ppp)
            {
                myfile << yprof[ppp] << "    " << uprof[ppp] << std::endl;
            }
        }
    }
    
    spade::convective::totani_lr tscheme(air);
    spade::convective::weno_3    wscheme(air);
    spade::viscous::visc_lr      visc_scheme(visc_law);
    
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
    const real_t time0    = 0.0;
    const real_t dx       = coords.ycoord.map(bounds.min(1)+grid.get_dx(1))-coords.ycoord.map(bounds.min(1));
    const real_t umax_ini = spade::algs::transform_reduce(prim, get_u, max_op);
    const real_t dt       = targ_cfl*dx/umax_ini;
    
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
        rhs = 0.0;
        grid.exchange_array(q);
        set_channel_noslip(q);
        const real_t alpha = 0.00005;
        spade::pde_algs::flux_div(q, rhs, wscheme);
        rhs *= (alpha)/(1.0-alpha);
        spade::pde_algs::flux_div(q, rhs, tscheme);
        rhs *= 1.0-alpha;
        spade::pde_algs::flux_div(q, rhs, visc_scheme);
        spade::pde_algs::source_term(rhs, [&]()->v5d{return v5d(0,0,force_term,0,0);});
    };
    
    spade::time_integration::rk2 time_int(prim, rhs0, rhs1, time0, dt, calc_rhs, ftrans, itrans);
    
    std::ofstream myfile("hist.dat");
    for (auto nti: range(0, nt_max))
    {
        int nt = nti;
        real_t umax     = spade::algs::transform_reduce(prim, get_u,       max_op);
        real_t ratio    = umax/(0.5*re_tau*u_tau + sqrt(air.gamma*air.R*t0));
        if (group.isroot())
        {
            const real_t cfl = umax*dt/dx;
            print(
                "nt:     ", spade::utils::pad_str(nt, 10),
                "cfl:    ", spade::utils::pad_str(cfl, 10),
                "u+a:    ", spade::utils::pad_str(umax, 10),
                "dx:     ", spade::utils::pad_str(dx, 10),
                "%:      ", spade::utils::pad_str(ratio, 10),
                "dt:     ", spade::utils::pad_str(dt, 10),
                "ftt:    ", spade::utils::pad_str(0.5*re_tau*u_tau*time_int.time()/delta, 10)
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
