#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
namespace cvdf::convective
{
    template <class T, class state_t>
    concept lr_conv_flux_function = fluid_state::is_state_type<state_t> && 
    requires(T t, state_t ql, state_t qr, const ctrs::array<typename state_t::value_type,3>& normal)
    {
        { t.calc_flux(ql, qr, normal) } -> cvdf::ctrs::basic_array;
    };
    
    template <fluid_state::state_dependent_gas gas_t> struct totani_lr
    {
        totani_lr(const gas_t& gas_in) {gas = &gas_in;}
        template <typename dtype> fluid_state::flux_t<dtype>
        calc_flux(
            const fluid_state::prim_t<dtype>& ql,
            const fluid_state::prim_t<dtype>& qr,
            const ctrs::array<dtype,3>& normal) const
        {
            fluid_state::flux_t<dtype> output;
            dtype un_l = normal[0]*ql.u()+normal[1]*ql.v()+normal[2]*ql.w();
            dtype un_r = normal[0]*qr.u()+normal[1]*qr.v()+normal[2]*qr.w();
            
            dtype rho_l = ql.p()/(gas->get_R(ql)*ql.T());
            dtype rho_r = qr.p()/(gas->get_R(qr)*qr.T());
            
            dtype e_l = ql.p()/(rho_l*(gas->get_gamma(ql)-1.0));
            dtype e_r = qr.p()/(rho_r*(gas->get_gamma(qr)-1.0));
            
            dtype c = 0.25*(rho_l+rho_r)*(un_l+un_r);
            dtype pm = 0.5*(ql.p() + qr.p());
            output.continuity() = c;
            output.energy()     = 0.5*c*(e_l + e_r + ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w()) + 0.5*(un_l*qr.p() + un_r*ql.p());
            output.x_momentum() = 0.5*c*(ql.u()+qr.u()) + normal[0]*pm;
            output.y_momentum() = 0.5*c*(ql.v()+qr.v()) + normal[1]*pm;
            output.z_momentum() = 0.5*c*(ql.w()+qr.w()) + normal[2]*pm;
            return output;
        }
        
        const gas_t* gas;
    };
}