#include <iostream>

template <typename thing_t> void print_type(const thing_t& t)
{
    //g++ only
    std::string pf(__PRETTY_FUNCTION__);
    std::size_t start = std::string("void print_type(const thing_t&) [with thing_t = ").length();
    std::size_t end = pf.length()-1;
    std::cout << pf.substr(start, end-start) << std::endl;
}

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
    real_t                             ffill = 0.0; //scalar eddy viscosity
    const real_t eps = 1.0;
    spade::state_sensor::ducros_t ducros(eps);
    using prim_t = decltype(pfill);
    
    std::vector<int> nx{8, 12, 24, 32, 40};
    
    std::vector<decltype(ffill)> er_2_all;
    std::vector<decltype(ffill)> er_i_all;
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
        
        enum syms {xs,ys,zs,uis,xis};
        
        symd::var_t<xs>  x;
        symd::var_t<ys>  y;
        symd::var_t<zs>  z;
        
        
        //Test functions
        const real_t alpha  = spade::consts::pi;
        auto p              = symd::cos(3.0*alpha*x)+symd::sin(2.0*alpha*y) + 10.0;
        auto T              = symd::one();
        auto u              = symd::cos(3.0*alpha*x)+symd::sin(2.0*alpha*y)*symd::sin(2.0*alpha*z);
        auto v              = symd::sin(3.0*alpha*x)-symd::cos(2.0*alpha*y)+symd::sin(2.0*alpha*z);
        auto w              = symd::cos(3.0*alpha*y)*symd::sin(2.0*alpha*x)+symd::sin(2.0*alpha*z);

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

        auto theta   = symd::ddx(u, x) + symd::ddx(v, y) + symd::ddx(w, z);
        auto theta2  = theta*theta;
        auto vorti0  = symd::ddx(w, y) - symd::ddx(v, z);
        auto vorti1  = symd::ddx(u, z) - symd::ddx(w, x);
        auto vorti2  = symd::ddx(v, x) - symd::ddx(u, y);
        auto vort2   = vorti0*vorti0 + vorti1*vorti1 + vorti2*vorti2;

        auto sens = theta2 / (theta2 + vort2 + eps);
        
        grid.exchange_array(prim);
        
        // shock sensor.
        auto num_alph = [&]()
        {
            using ss_t = decltype(ducros);
            auto kernel = spade::omni::make_kernel([](const ss_t& dc, const auto& info){return dc.get_sensor(info);}, ducros);
            spade::algs::transform_to(prim, rhs, kernel);
            spade::io::output_vtk("output", "num", rhs);
        };
        auto ana_alph = [&]()
        {
            spade::algs::fill_array(rhs_ana, [=](const grid_type::coord_point_type& xg)
            {
                return sens(x=xg[0], y=xg[1], z=xg[2]);
            });
            spade::io::output_vtk("output", "ana", rhs_ana);
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
            real_t l2_error;
            real_t li_error;
            spade::reduce_ops::reduce_max<real_t> max_op;
            spade::reduce_ops::reduce_sum<real_t> sum_op;
            li_error = spade::algs::transform_reduce(rhs, [=](const real_t& rhsi) {return spade::utils::abs(rhsi);}, max_op);
            l2_error = std::sqrt(dV*spade::algs::transform_reduce(rhs, [=](const real_t& rhsi) {return rhsi*rhsi;}, sum_op));
            return std::make_pair(l2_error, li_error);
        };
        
        auto [er_2, er_i] = run_test(num_alph, ana_alph);
        
        dx_all.push_back(std::pow(dV, 1.0/3.0));
        er_2_all.push_back(er_2);
        er_i_all.push_back(er_i);
    }

    std::vector<decltype(ffill)> er_2_ord;
    std::vector<decltype(ffill)> er_i_ord;
    
    for (int l = 0; l < nx.size()-1; ++l)
    {
        auto comp_ord = [](auto y1, auto y0, auto x1, auto x0)
        {
            auto r = y1;
            r = (std::log(y1) - std::log(y0)) / (std::log(x1) - std::log(x0));
            return r;
        };
        auto l_er_2 = comp_ord(er_2_all[l+1], er_2_all[l], dx_all[l+1], dx_all[l]);
        auto l_er_i = comp_ord(er_i_all[l+1], er_i_all[l], dx_all[l+1], dx_all[l]);
        
        er_2_ord.push_back(l_er_2);
        er_i_ord.push_back(l_er_i);
    }
    
    
    if (group.isroot())
    {
        print("=========== Eddy visc. ===========");
        print(">> L2:");
        for (auto oo: er_2_ord) print("  ", oo);
        print(">> LInf:");
        for (auto oo: er_i_ord) print("  ", oo);
    }
    return 0;
}
