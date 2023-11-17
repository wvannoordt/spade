#pragma once

#include "dispatch/shared_vector.h"
#include "sym/symbol.h"

namespace spade::ode
{
    template <sym::is_symbol sym_t, typename data_t>
    struct constant_mesh_t
    {
        using symbol_type = sym_t;
        using vector_type = device::shared_vector<data_t>;
        vector_type xs;
        std::size_t count;
        
        constant_mesh_t(const sym_t&, const vector_type& xs_in, const std::size_t count_in)
        : xs{xs_in}, count{count_in} {}
        
        std::size_t num_instances() const { return count; }
        std::size_t total_elements() const { return num_instances()*xs.size(); }
        
        data_t get_x(const int i, const int) const
        {
            return xs[i];
        }
    };
}