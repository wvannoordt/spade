#pragma once

#include "ode/expression.h"

namespace spade::ode
{
    template <typename osystem_t, typename...exprs_t>
    // requires(something complicated) 
    static void solve_bvp(osystem_t& sys, const expression_list_t<exprs_t...> expressions)
    {
        using expr_list_type = expression_list_t<exprs_t...>;
        using vars_t         = typename expr_list_type::syms_type;
        const std::size_t mesh_size = sys.mesh.mesh_size();
        const std::size_t instances = sys.mesh.num_instances();
        
        dispatch::ranges::bilinear_range_t rng(0UL, mesh_size, 0UL, instances, sys.device());
        using index_type = decltype(rng)::index_type;
        
        auto img = sys.image();
        
        //Initialize according to the linear profiles
        auto loop = [=] _sp_hybrid (const index_type& i) mutable
        {
            //Instance index
            const auto inst = i[1];
            
            //Wall index
            const auto iw   = i[0];
            auto buf = img.get_buffer(inst);
            algs::static_for<0, vars_t::size()>([&](const auto& jj)
            {
                constexpr static int j = jj.value;
                using sym_t      = typename vars_t::elem<j>;
                const auto expr = get_expression(sym_t(), expressions);
                
                // If expr is an ODE, this gives the boundary conditions.
                // If not an ODE, gives the kind of expression that expr is (explicit, algebraic, etc)
                const auto spec      = expr.specifier;
                const auto mesh_sym  = sys.mesh_symbol();
                if constexpr (is_ode<decltype(spec)>)
                {
                    const auto m_b  = detail::lin_coeffs(spec, buf[mesh_sym][0], buf[mesh_sym].back());
                    buf[sym_t()][iw] = m_b[0] + buf[mesh_sym][iw]*m_b[1];
                }
            });
        };
        
        dispatch::execute(rng, loop);
    }
}