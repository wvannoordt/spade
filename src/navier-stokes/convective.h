#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
namespace cvdf::convective
{
    template <class T, const std::size_t dim, class state_t, class gas_t>
    concept lr_conv_flux_function = fluid_state::is_state_type<state_t> && fluid_state::state_dependent_gas<T> && 
    requires(T t, state_t ql, state_t qr, gas_t gas, ctrs::array<typename state_t::value_type,3> normal)
    {
        { t.calc_flux(ql, qr, normal, gas) } -> cvdf::ctrs::vec_nd<dim, typename state_t::value_type>;
    };
    
    struct totani_lr
    {
        template <typename dtype, fluid_state::state_dependent_gas gas_t> ctrs::array<dtype,5>
        calc_flux(
            const fluid_state::prim_t<dtype>& ql,
            const fluid_state::prim_t<dtype>& qr,
            const ctrs::array<dtype,3>& normal,
            const gas_t& gas)
        {
            return ctrs::array<dtype,5>(0,0,0,0,0);
        }
    };
}