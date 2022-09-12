#include <chrono>
#include "spade.h"

typedef double real_t;
typedef spade::ctrs::array<real_t, 3> v3d;
typedef spade::ctrs::array<int,    3> v3i;
typedef spade::ctrs::array<int,    4> v4i;
typedef spade::ctrs::array<spade::grid::cell_t<int>, 4> v4c;
typedef spade::fluid_state::prim_t<real_t> prim_t;
typedef spade::fluid_state::flux_t<real_t> flux_t;
typedef spade::fluid_state::cons_t<real_t> cons_t;

#include "calc_u_bulk.h"
#include "calc_boundary_flux.h"

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
                        spade::grid::cell_idx_t i_d(ii[0], j-(nnn+0)*nvec_out[1], ii[1], lb);
                        spade::grid::cell_idx_t i_g(ii[0], j+(nnn+1)*nvec_out[1], ii[1], lb);
                        prim_t q_d, q_g;
                        for (auto n: range(0,5)) q_d[n] = prims(n, i_d[0], i_d[1], i_d[2], i_d[3]);
                        const auto x_g = grid.get_comp_coords(i_g);
                        const auto x_d = grid.get_comp_coords(i_d);
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
            real_t d2p = prims(2+idir, i[0],   i[1],   i[2]+1, i[3]);
            real_t d2m = prims(2+idir, i[0],   i[1],   i[2]-1, i[3]);
            lap += (d0p - 2*d00 + d0m)/(dx*dx);
            lap += (d1p - 2*d00 + d1m)/(dy*dy);
            lap += (d2p - 2*d00 + d2m)/(dz*dz);
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
    
    bool viz_only = false;
    std::string viz_filename = "";
    if (args.has_arg("-viz"))
    {
        viz_filename = args["-viz"];
        if (group.isroot()) print("Converting", viz_filename);
        viz_only = true;
        if (!std::filesystem::exists(viz_filename))
        {
            if (group.isroot()) print("Cannot find viz file", viz_filename);
            abort();
        }
    }
    
    spade::ctrs::array<int, 3> num_blocks(8, 8, 8);
    spade::ctrs::array<int, 3> cells_in_block(48, 48, 48);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    //spade::ctrs::array<int, 3> exchange_cells(8, 8, 8);
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
    const int    nt_max   = 300001;
    const int    nt_skip  = 50000000;
    const int    checkpoint_skip  = 5000;
    
    spade::coords::identity<real_t> coords;
    
    std::filesystem::path out_path("checkpoint");
    if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
    
    
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    
    prim_t fill1 = 0.0;
    flux_t fill2 = 0.0;
    
    spade::grid::grid_array prim (grid, fill1);
    
    if (viz_only)
    {
        spade::io::binary_read(viz_filename, prim);
        std::string fname = spade::io::output_vtk("output", "viz", prim);
        if (group.isroot()) print("outputting", fname);
        return 0;
    }
    
    spade::grid::grid_array rhs (grid, fill2);
    
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
    
    const int nidx = 8;
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
        if (group.isroot()) print("reading...");
        spade::io::binary_read(init_filename, prim);
        if (group.isroot()) print("Init done.");
        grid.exchange_array(prim);
        set_channel_noslip(prim);
    }
    
    spade::convective::totani_lr tscheme(air);
    spade::convective::weno_3    wscheme(air);
    spade::viscous::visc_lr  visc_scheme(visc_law);
    
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
    
    
    
    struct trans_t
    {
        using gas_t = spade::fluid_state::perfect_gas_t<real_t>;
        const gas_t* gas;
        
        struct p2c_t
        {
            const gas_t* gas;
            typedef prim_t arg_type;
            p2c_t(const gas_t& gas_in) {gas = &gas_in;}
            cons_t operator () (const prim_t& q) const
            {
                cons_t w;
                spade::fluid_state::convert_state(q, w, *gas);
                return w;
            };
        };
        struct c2p_t
        {
            const gas_t* gas;
            typedef cons_t arg_type;
            c2p_t(const gas_t& gas_in) {gas = &gas_in;}
            prim_t operator () (const cons_t& w) const
            {
                prim_t q;
                spade::fluid_state::convert_state(w, q, *gas);
                return q;
            }
        };
        
        trans_t(const gas_t& gas_in) { gas = &gas_in; }
        
        void transform_forward (decltype(prim)& q) const { spade::algs::transform_inplace(q, p2c_t(*gas)); }
        void transform_inverse (decltype(prim)& q) const { spade::algs::transform_inplace(q, c2p_t(*gas)); }
    } trans(air);
    
    
    auto calc_rhs = [&](auto& rhs, auto& q, const auto& t) -> void
    {
        rhs = 0.0;
        grid.exchange_array(q);
        set_channel_noslip(q);
        spade::pde_algs::flux_div(q, rhs, tscheme);
        spade::pde_algs::flux_div(prim, rhs, visc_scheme);
        spade::algs::transform_inplace(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_ar) -> spade::ctrs::array<real_t, 5> 
        {
            spade::ctrs::array<real_t, 5> rhs_new = rhs_ar;
            rhs_new[2] += force_term;
            return rhs_new;
        });
    };
    
    spade::time_integration::rk2 time_int(prim, rhs, time0, dt, calc_rhs, trans);
    
    std::ofstream myfile("hist.dat");
    for (auto nti: range(0, nt_max))
    {
        int nt = nti;
        const real_t umax   = spade::algs::transform_reduce(prim, get_u, max_op);
        real_t ub, rhob;
        calc_u_bulk(prim, air, ub, rhob);
        const real_t area = bounds.size(0)*bounds.size(2);
        auto conv2 = proto::get_domain_boundary_flux(prim, visc_scheme, 2);
        auto conv3 = proto::get_domain_boundary_flux(prim, visc_scheme, 3);
        conv2 /= area;
        conv3 /= area;
        const real_t tau = 0.5*(spade::utils::abs(conv2.x_momentum()) + spade::utils::abs(conv3.x_momentum()));
        
        if (group.isroot())
        {
            const real_t cfl = umax*dt/dx;
            const int pn = 10;
            print(
                "nt: ", spade::utils::pad_str(nt, pn),
                "cfl:", spade::utils::pad_str(cfl, pn),
                "u+a:", spade::utils::pad_str(umax, pn),
                "ub: ", spade::utils::pad_str(ub, pn),
                "rb: ", spade::utils::pad_str(rhob, pn),
                "tau:", spade::utils::pad_str(tau, pn),
                "dx: ", spade::utils::pad_str(dx, pn),
                "dt: ", spade::utils::pad_str(dt, pn),
                "ftt:", spade::utils::pad_str(20.0*u_tau*time_int.time()/delta, pn)
            );
            myfile << nt << " " << cfl << " " << umax << " " << ub << " " << rhob << " " << tau << " " << dx << " " << dt << std::endl;
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
    	auto start = std::chrono::steady_clock::now();
        time_int.advance();
    	auto end = std::chrono::steady_clock::now();
    	if (group.isroot()) print("Elapsed:", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), "ms");
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
