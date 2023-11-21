#pragma once

#include "ode/expression.h"

namespace spade::ode
{
    template <typename osystem_t, typename...exprs_t>
    // requires(something complicated) 
    static void solve_bvp(const osystem_t& sys, const expression_list_t<exprs_t...> expressions)
    {
        using expr_list_type = expression_list_t<exprs_t...>;
        const std::size_t mesh_size = sys.mesh.mesh_size();
        const std::size_t instances = sys.mesh.num_instances();
        
        dispatch::ranges::bilinear_range_t rng(0UL, mesh_size, 0UL, instances, sys.device());
        using index_type = decltype(rng)::index_type;
        
        const auto img = sys.image();
        
        auto loop = [=] _sp_hybrid (const index_type& i) mutable
        {
            const auto inst = i[1];
            const auto buf  = img.get_buffer(inst);
        };
        
        dispatch::execute(rng, loop);
    }
}