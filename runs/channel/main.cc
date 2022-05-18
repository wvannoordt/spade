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
    const real_t t_wall = 325.0;
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
                    v4c i_d(ii[0], j,             ii[1], lb[0]);
                    v4c i_g(ii[0], j+nvec_out[1], ii[1], lb[0]);
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
    const double duhat = 0.5;
    auto ini = [&](const cvdf::ctrs::array<real_t, 3> x, const int& i, const int& j, const int& k, const int& lb) -> prim_t
    {
        prim_t output;
        output.p() = p0;
        output.T() = t0;
        output.u() = 1.5*(1.0 - x[1]*x[1]/(delta*delta))*u0;
        output.v() = 0.0;
        output.w() = 0.0;
        if ((((i/6) + (j/6))%2==0) && (((i/3) + (j/3))%2==0))
        {
            output.u() += duhat*u0;
            output.v() += duhat*u0;
            output.w() += duhat*u0;
        }
        return output;
    };
    
    cvdf::algs::fill_array(prim, ini);
    cvdf::output::output_vtk("output", "ini", grid, prim);
    cvdf::convective::totani_lr tscheme(air);
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
    
    const int    nt_max = 100000;
    const real_t dt     = 1e-2;
    for (auto nti: range(0, nt_max))
    {
        cvdf::algs::transform_inplace(prim, c2p);
        int nt = nti[0];
        grid.exchange_array(prim);
        set_channel_noslip(prim);
        cvdf::flux_algs::flux_lr_diff(prim, rhs, tscheme);
        // cvdf::flux_algs::flux_lr_diff(prim, rhs, visc_scheme);
        cvdf::algs::transform_inplace(prim, p2c);
        rhs  *= dt;
        prim += rhs;
        real_t umax   = cvdf::algs::transform_reduce(prim, get_u, max_op);
        real_t rhsmax = cvdf::algs::transform_reduce(rhs,  [](const cvdf::ctrs::array<real_t, 5>& rhs) -> real_t {return cvdf::utils::max(rhs[0], rhs[1], rhs[2], rhs[3], rhs[4]);}, max_op);
        if (group.isroot())
        {
            print(nt, umax, rhsmax);
        }
        if (nt%300 == 0)
        {
            cvdf::output::output_vtk("output", "last", grid, prim);
        }
    }
    
    
    return 0;
}