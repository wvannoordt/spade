#pragma once

#include "core/converge_crit.h"
#include "ode/expression.h"

namespace spade::ode
{
    template <typename osystem_t, typename...exprs_t>
    // requires(something complicated) 
    static void solve_bvp(
        osystem_t& sys,
        const expression_list_t<exprs_t...> expressions,
        const utils::converge_crit_t<typename osystem_t::value_type>& params)
    {
        using expr_list_type = expression_list_t<exprs_t...>;
        using vars_t         = typename expr_list_type::syms_type;
        const std::size_t mesh_size = sys.mesh.mesh_size();
        const std::size_t instances = sys.mesh.num_instances();
        
        //Really should put together something better than this...
        dispatch::ranges::bilinear_range_t rng(0UL, mesh_size, 0UL, instances, sys.device());
        using index_type = decltype(rng)::index_type;
        
        auto img = sys.image();
        
        //Initialize according to the linear profiles
        auto init = [=] _sp_hybrid (const index_type& i) mutable
        {
            //Instance index
            const std::size_t inst = std::size_t(i[1]);
            
            //Wall index
            const int iw   = int(i[0]);
            auto buf = img.get_buffer(inst);
            
            //Initialize the each variable (don't do this if we want to keep initial guess)
            algs::static_for<0, vars_t::size()>([&](const auto& jj)
            {
                constexpr static int j = jj.value;
                using sym_t       = typename vars_t::elem<j>;
                const auto& expr  = get_expression(sym_t(), expressions);
                
                // If expr is an ODE, this gives the boundary conditions.
                // If not an ODE, gives the kind of expression that expr is (explicit, algebraic, etc)
                const auto spec      = expr.specifier;
                const auto mesh_sym  = sys.mesh_symbol();
                
                //Initialize the solution with linear initial guess
                if constexpr (is_ode<decltype(spec)>)
                {
                    const auto m_b   = detail::lin_coeffs(spec, buf[mesh_sym][0], buf[mesh_sym].back());
                    buf[sym_t()][iw] = m_b[0] + buf[mesh_sym][iw]*m_b[1];
                }
                else
                {
                    buf[sym_t()][iw] = expr.func(inst, iw, buf);
                }
            });
        };
        
        using real_t = typename osystem_t::value_type;
        const auto& device = sys.device();
        
        constexpr int num_vars = osystem_t::variable_list::size();
        ctrs::array<real_t, num_vars> eps, eps_i;
        
        // dispatch::shared_t shdata(eps_i, eps);
        
        // auto inner_range = dispatch::ranges::make_range({0UL, mesh_size}, device);
        // auto outer_range = dispatch::ranges::make_range({0UL, instances}, device);
        // dispatch::threads_t g_threads(inner_range);
        
        // auto solve = [=] _sp_hybrid (const std::size_t& i_inst, const dispatch::threads_t& thrds) mutable
        // {
        //     auto& eps_ini = data[0_c];
        //     auto& eps_cur = data[1_c];
            
        //     auto buf = img.get_buffer(inst);
            
        //     if (thrds.isroot()) eps_cur = real_t(0.0);
        //     if (thrds.isroot()) eps_ini = real_t(0.0);
            
        //     thrds.do([&](const int iw)
        //     {
                
        //     });
            
            
            
        //     thrds.sync();
            
        //     int its = 0;
        //     while(its++ < params.max_its())
        //     {
        //         // Newton Method:
        //         // x_{n+1} = x_{n} - J^{-1}(x_{n})*f(x_{n})
                
        //         //First, update all of the algebraic relations:
                
                
        //         //Need some kind of check
        //     }
        // };
        
        // dispatch::kernel_t k_solve(inner, shdata);
        
        dispatch::execute(rng, init);
        // dispatch::execute(outer_range, g_threads, k_solve);
    }
}