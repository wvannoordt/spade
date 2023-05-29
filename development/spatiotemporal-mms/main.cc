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

        const real_t c0 = -1.0/12.0;
        const real_t c1 =  7.0/12.0;
        const real_t c2 =  7.0/12.0;
        const real_t c3 = -1.0/12.0;
        return u*(c0*p0 + c1*p1 + c2*p2 + c3*p3);
    }
};

template <typename in_t> auto f_mod(const in_t& inp, const in_t& modu)
{                                                                                                            
    if (inp>0) return fmod(inp, modu);
    return fmod(inp + ((int)(-inp/modu) + 1)*modu, modu);
}

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);

    std::vector<int> nx{8, 16, 22, 30, 48};
    std::vector<real_t> dx_v;
    std::vector<real_t> er_v;
    for (int lev = 0; lev < nx.size(); ++lev)
    {
        if (group.isroot()) print("nx =", nx[lev]);
        spade::ctrs::array<int, 2> num_blocks(2, 2);
        spade::ctrs::array<int, 2> cells_in_block(nx[lev], nx[lev]);
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
        spade::grid::grid_array phia(grid, fill1);
        spade::grid::grid_array rhs (grid, fill1);

        using point_type = decltype(grid)::coord_point_type;

        auto ini = [&](const point_type& x)
        {
            const real_t r = std::sqrt(x[1]*x[1]+x[0]*x[0]);
            return std::exp(-15.0*r*r);
        };
        
        spade::algs::fill_array(phi, ini);
        grid.exchange_array(phi);
        
        spade::ctrs::array<real_t, 2> vel(0.3, 1.6);
        const real_t targ_cfl = 0.5;
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
        spade::time_integration::time_axis_t       axis    (time0, dt);
        spade::time_integration::rk4_t             alg;
        spade::time_integration::integrator_data_t q       (phi, rhs, alg);
        spade::time_integration::integrator_t      time_int(axis, alg, q, calc_rhs, boundary_cond);

        const real_t tmax = 2.0*bounds.size(0)/umax_ini;

        spade::timing::mtimer_t tmr("advance");
        const int nt_max  = 10000;
        const int nt_skip = nt_max/100;
        int nt = 0;
        spade::reduce_ops::reduce_max<real_t> time_max;
        spade::reduce_ops::reduce_max<real_t> max_op;
        const auto linf     = [](const real_t& val){ return spade::utils::abs(val); };
        time_max.init(0.0);

        while (time_int.time() < tmax)
        {
            const auto& sol = time_int.solution();
            const real_t tm = time_int.time();
            spade::algs::fill_array(phia, [&](const point_type& x)
            {
                real_t xi  = fmod(x[0] - vel[0]*tm + bounds.min(0), bounds.size(0)) - bounds.min(0);
                real_t eta = fmod(x[1] - vel[1]*tm + bounds.min(1), bounds.size(1)) - bounds.min(1);
                return ini(point_type(xi, eta));
            });

            phia -= sol;

            const real_t er_max = spade::algs::transform_reduce(phia, linf, max_op);
            time_max.reduce_elem(er_max);

            tmr.start("advance");
            time_int.advance();
            tmr.stop("advance");
            
            ++nt;
        }
        dx_v.push_back(grid.get_dx(0));
        er_v.push_back(time_max.get_value());
    }

    std::vector<real_t> ords;
    for (int p = 0; p < dx_v.size()-1; ++p)
    {
        ords.push_back((std::log(er_v[p+1]) - std::log(er_v[p])) / (std::log(dx_v[p+1]) - std::log(dx_v[p])));
    }

    for (auto o: ords)
    {
        if (group.isroot()) print("spatiotemporal order:", o);
    }
    
    return 0;
}
