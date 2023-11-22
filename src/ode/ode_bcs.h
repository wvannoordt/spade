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
    
    namespace detail
    {
        template <typename data_t, typename dx_t>
        _sp_hybrid static constexpr ctrs::array<data_t, 2>
        lin_coeffs(const ode_bc_t<neumann_t<data_t>, dirichlet_t<data_t>>& bcs, const dx_t& x0, const dx_t& x1)
        {
            return {bcs.right.value - bcs.left.value*x1, bcs.left.value};
        }
        
        template <typename data_t, typename dx_t>
        _sp_hybrid static constexpr ctrs::array<data_t, 2>
        lin_coeffs(const ode_bc_t<dirichlet_t<data_t>, neumann_t<data_t>>& bcs, const dx_t& x0, const dx_t& x1)
        {
            return {bcs.left.value - bcs.right.value*x0, bcs.right.value};
        }
        
        template <typename data_t, typename dx_t>
        _sp_hybrid static constexpr ctrs::array<data_t, 2>
        lin_coeffs(const ode_bc_t<dirichlet_t<data_t>, dirichlet_t<data_t>>& bcs, const dx_t& x0, const dx_t& x1)
        {
            const auto dxinv = data_t(1.0)/(x1-x0);
            return {dxinv*x1*bcs.left.value - dxinv*x0*bcs.right.value, (bcs.right.value - bcs.left.value)*dxinv};
        }
    }
}