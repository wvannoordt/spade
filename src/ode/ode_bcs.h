#pragma once

#include <concepts>

namespace spade::ode
{
    // "explicit" is a reserved keyword :(
    constexpr static struct explicit_t{} xplicit;
    
    template <typename T> concept is_ode = T::is_ode_bc_v;
    
    template <typename left_t, typename right_t>
    struct ode_bc_t
    {
        static constexpr bool is_ode_bc_v = true;
        left_t  left;
        right_t right;
    };
    
    static struct dirichlet_t { } dirichlet;
    static struct zerograd_t  { } zerograd;
    
    template <typename left_t, typename right_t>
    static auto make_bcs(const left_t& left, const right_t& right)
    {
        return ode_bc_t{left, right};
    }
}