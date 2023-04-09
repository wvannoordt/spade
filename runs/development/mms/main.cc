#include "symd.h"
#include "spade.h"

int main(int argc, char** argv)
{
    using real_t = double;
    spade::parallel::mpi_t group(&argc, &argv);
    
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) = 0.0;
    bounds.max(0) = 2.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 2.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 2.0;
    spade::coords::identity<real_t> coords;
    spade::fluid_state::prim_t<real_t> pfill = 0.0;
    spade::fluid_state::flux_t<real_t> ffill = 0.0;
    
    std::vector<int> nx{8, 12, 24, 32, 40};
    
    std::vector<decltype(ffill)> c_er_2_all;
    std::vector<decltype(ffill)> c_er_i_all;
    std::vector<decltype(ffill)> v_er_2_all;
    std::vector<decltype(ffill)> v_er_i_all;
    std::vector<real_t> dx_all;
    for (int l = 0; l < nx.size(); ++l)
    {
        
        spade::ctrs::array<int, 3> cells_in_block(nx[l], nx[l], nx[l]);
        spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
        const auto dV = grid.get_dx(0)*grid.get_dx(1)*grid.get_dx(2);

        spade::grid::grid_array prim   (grid, pfill);
        spade::grid::grid_array rhs    (grid, ffill);
        spade::grid::grid_array rhs_ana(grid, ffill);
        using grid_type = decltype(grid);
        
        enum syms {xs,ys,zs};
        
        symd::var_t<xs> x;
        symd::var_t<ys> y;
        symd::var_t<zs> z;
        
        const real_t alpha   = spade::consts::pi;
        const real_t rgas    = 287.15;
        const real_t gamma   = 1.4;
        const real_t mu      = 1.0;
        const real_t prandtl = 0.72;
        const real_t cond    = (rgas*gamma/(gamma-1.0))*mu/prandtl;
        
        //Test functions
        auto p   = 5.0+2.0*symd::cos(3.0*alpha*x)*symd::sin(2.0*alpha*y)+symd::sin(4.0*alpha*z);
        auto T   = 10.0+2.0*symd::cos(2.0*alpha*x)*symd::sin(3.0*alpha*y)+symd::sin(4.0*alpha*z);
        auto u   = symd::sin(3.0*alpha*x)*symd::cos(2.0*alpha*y)*symd::cos(2.0*alpha*z);
        auto v   = symd::cos(3.0*alpha*x)*symd::cos(2.0*alpha*y)*symd::cos(3.0*alpha*z);
        auto w   = symd::sin(3.0*alpha*x)*symd::sin(2.0*alpha*y)*symd::cos(4.0*alpha*z);
        
        auto rho = p/(rgas*T);
        auto Et  = (rgas*T/(gamma-1.0)) + 0.5*(u*u + v*v + w*w);
        
        auto div = symd::ddx(u, x) + symd::ddx(v, y) + symd::ddx(w, z);
        
        auto txx = mu*(symd::ddx(u, x) + symd::ddx(u, x)) - 0.66666666667*mu*div;
        auto txy = mu*(symd::ddx(u, y) + symd::ddx(v, x));
        auto txz = mu*(symd::ddx(u, z) + symd::ddx(w, x));
        
        auto tyx = mu*(symd::ddx(v, x) + symd::ddx(u, y));
        auto tyy = mu*(symd::ddx(v, y) + symd::ddx(v, y)) - 0.66666666667*mu*div;
        auto tyz = mu*(symd::ddx(v, z) + symd::ddx(w, y));
        
        auto tzx = mu*(symd::ddx(w, x) + symd::ddx(u, z));
        auto tzy = mu*(symd::ddx(w, y) + symd::ddx(v, z));
        auto tzz = mu*(symd::ddx(w, z) + symd::ddx(w, z)) - 0.66666666667*mu*div;
        
        auto qx = -cond*symd::ddx(T, x);
        auto qy = -cond*symd::ddx(T, y);
        auto qz = -cond*symd::ddx(T, z);

        spade::algs::fill_array(prim, [=](const grid_type::coord_point_type& xg)
        {
            spade::fluid_state::prim_t<real_t> output;
            output.p() = p(x=xg[0], y=xg[1], z=xg[2]);
            output.T() = T(x=xg[0], y=xg[1], z=xg[2]);
            output.u() = u(x=xg[0], y=xg[1], z=xg[2]);
            output.v() = v(x=xg[0], y=xg[1], z=xg[2]);
            output.w() = w(x=xg[0], y=xg[1], z=xg[2]);
            return output;
        });
        
        grid.exchange_array(prim);
        
        auto cont_i  = -1*(symd::ddx(rho*u,          x) + symd::ddx(rho*v,          y) + symd::ddx(rho*w,          z));
        auto mom_x_i = -1*(symd::ddx(rho*u*u + p,    x) + symd::ddx(rho*v*u,        y) + symd::ddx(rho*w*u,        z));
        auto mom_y_i = -1*(symd::ddx(rho*u*v,        x) + symd::ddx(rho*v*v + p,    y) + symd::ddx(rho*w*v,        z));
        auto mom_z_i = -1*(symd::ddx(rho*u*w,        x) + symd::ddx(rho*v*w,        y) + symd::ddx(rho*w*w + p,    z));
        auto eng_i   = -1*(symd::ddx(rho*u*Et + u*p, x) + symd::ddx(rho*v*Et + v*p, y) + symd::ddx(rho*w*Et + w*p, z));
        
        //no cont_v
        auto mom_x_v = symd::ddx(txx, x) + symd::ddx(txy, y) + symd::ddx(txz, z);
        auto mom_y_v = symd::ddx(tyx, x) + symd::ddx(tyy, y) + symd::ddx(tyz, z);
        auto mom_z_v = symd::ddx(tzx, x) + symd::ddx(tzy, y) + symd::ddx(tzz, z);
        auto eng_v   = symd::ddx(u*txx+v*txy+w*txz-qx, x) + symd::ddx(u*tyx+v*tyy+w*tyz-qy, y) + symd::ddx(u*tzx+v*tzy+w*tzz-qz, z);
        
        spade::fluid_state::ideal_gas_t<real_t> air;
        air.R = rgas;
        air.gamma = gamma;
        spade::viscous_laws::constant_viscosity_t visc_law(mu, prandtl, air);

        spade::convective::totani_lr tscheme(air);
        spade::viscous::visc_lr  visc_scheme(visc_law);
        
        // Convective
        auto num_conv = [&](){spade::pde_algs::flux_div(prim, rhs, tscheme);};
        auto ana_conv = [&]()
        {
            spade::algs::fill_array(rhs_ana, [=](const grid_type::coord_point_type& xg)
            {
                spade::fluid_state::flux_t<real_t> output;
                output.continuity() = cont_i (x=xg[0], y=xg[1], z=xg[2]);
                output.energy()     = eng_i  (x=xg[0], y=xg[1], z=xg[2]);
                output.x_momentum() = mom_x_i(x=xg[0], y=xg[1], z=xg[2]);
                output.y_momentum() = mom_y_i(x=xg[0], y=xg[1], z=xg[2]);
                output.z_momentum() = mom_z_i(x=xg[0], y=xg[1], z=xg[2]);
                return output;
            });
        };
        
        // Viscous
        auto num_visc = [&](){spade::pde_algs::flux_div(prim, rhs, visc_scheme);};
        auto ana_visc = [&]()
        {
            spade::algs::fill_array(rhs_ana, [=](const grid_type::coord_point_type& xg)
            {
                spade::fluid_state::flux_t<real_t> output;
                output.continuity() = 0.0;
                output.energy()     = eng_v  (x=xg[0], y=xg[1], z=xg[2]);
                output.x_momentum() = mom_x_v(x=xg[0], y=xg[1], z=xg[2]);
                output.y_momentum() = mom_y_v(x=xg[0], y=xg[1], z=xg[2]);
                output.z_momentum() = mom_z_v(x=xg[0], y=xg[1], z=xg[2]);
                return output;
            });
        };
        
        auto run_test = [&](const auto& num_func, const auto& ana_func)
        {
            rhs = 0.0;
            rhs_ana = 0.0;
            
            spade::timing::mtimer_t tmr("num", "ana");
            
            tmr.start("num");
            num_func();
            tmr.stop("num");
            
            tmr.start("ana");
            ana_func();
            tmr.stop("ana");
            
            if (group.isroot()) print(tmr);
            rhs -= rhs_ana;
            grid.exchange_array(rhs);
            spade::fluid_state::flux_t<real_t> l2_error;
            spade::fluid_state::flux_t<real_t> li_error;
            spade::reduce_ops::reduce_max<real_t> max_op;
            spade::reduce_ops::reduce_sum<real_t> sum_op;
            for (int i = 0; i < l2_error.size(); ++i)
            {
                li_error[i] = spade::algs::transform_reduce(rhs, [=](const spade::fluid_state::flux_t<real_t>& rhsi) {return spade::utils::abs(rhsi[i]);}, max_op);
                l2_error[i] = std::sqrt(dV*spade::algs::transform_reduce(rhs, [=](const spade::fluid_state::flux_t<real_t>& rhsi) {return rhsi[i]*rhsi[i];}, sum_op));
            }
            return std::make_pair(l2_error, li_error);
        };
        
        auto [c_er_2, c_er_i] = run_test(num_conv, ana_conv);
        auto [v_er_2, v_er_i] = run_test(num_visc, ana_visc);
        
        dx_all.push_back(std::pow(dV, 1.0/3.0));
        c_er_2_all.push_back(c_er_2);
        c_er_i_all.push_back(c_er_i);
        v_er_2_all.push_back(v_er_2);
        v_er_i_all.push_back(v_er_i);
    }
    std::vector<decltype(ffill)> c_er_2_ord;
    std::vector<decltype(ffill)> c_er_i_ord;
    std::vector<decltype(ffill)> v_er_2_ord;
    std::vector<decltype(ffill)> v_er_i_ord;
    
    for (int l = 0; l < nx.size()-1; ++l)
    {
        auto comp_ord = [](auto y1, auto y0, auto x1, auto x0)
        {
            auto r = y1;
            for (int p = 0; p < r.size(); ++p)
            {
                r[p] = (std::log(y1[p]) - std::log(y0[p])) / (std::log(x1) - std::log(x0));
            }
            return r;
        };
        auto l_c_er_2 = comp_ord(c_er_2_all[l+1], c_er_2_all[l], dx_all[l+1], dx_all[l]);
        auto l_c_er_i = comp_ord(c_er_i_all[l+1], c_er_i_all[l], dx_all[l+1], dx_all[l]);
        auto l_v_er_2 = comp_ord(v_er_2_all[l+1], v_er_2_all[l], dx_all[l+1], dx_all[l]);
        auto l_v_er_i = comp_ord(v_er_i_all[l+1], v_er_i_all[l], dx_all[l+1], dx_all[l]);
        
        //No mass diffusion
        l_v_er_2.continuity() = 0.0;
        l_v_er_i.continuity() = 0.0;
        
        c_er_2_ord.push_back(l_c_er_2);
        c_er_i_ord.push_back(l_c_er_i);
        v_er_2_ord.push_back(l_v_er_2);
        v_er_i_ord.push_back(l_v_er_i);
    }
    
    
    if (group.isroot())
    {
        print("=========== Convective ===========");
        print(">> L2:");
        for (auto oo: c_er_2_ord) print("  ", oo);
        print(">> LInf:");
        for (auto oo: c_er_i_ord) print("  ", oo);
        
        
        print("===========  Viscous   ===========");
        print(">> L2:");
        for (auto oo: v_er_2_ord) print("  ", oo);
        print(">> LInf:");
        for (auto oo: v_er_i_ord) print("  ", oo);
    }
    return 0;
}
