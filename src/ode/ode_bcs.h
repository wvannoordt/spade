#pragma once

namespace spade::ode
{
    // "explicit" is a reserved keyword :(
    constexpr static struct explicit_t{} xplicit;
    
    template <typename left_t, typename right_t>
    struct ode_bc_t
    {
        left_t  left;
        right_t right;
    };
    
    template <typename val_t>
    struct neumann_t   { val_t value; };
    
    template <typename val_t>
    static neumann_t<val_t> neumann(const val_t& val)
    {
        return neumann_t{val};
    }
    
    template <typename val_t>
    struct dirichlet_t { val_t value; };
    
    template <typename val_t>
    static dirichlet_t<val_t> dirichlet(const val_t& val)
    {
        return dirichlet_t{val};
    }
    
    template <typename left_t, typename right_t>
    static auto make_bcs(const left_t& left, const right_t& right)
    {
        return ode_bc_t{left, right};
    }
}