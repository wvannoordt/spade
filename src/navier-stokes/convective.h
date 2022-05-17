#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/flux_base.h"
#include "navier-stokes/flux_input.h"

namespace cvdf::convective
{
    template <class T, class state_t>
    concept lr_conv_flux_function = fluid_state::is_state_type<state_t> && 
    requires(
        T t,
        state_t ql,
        state_t qr,
        const ctrs::array<typename state_t::value_type,3>& normal_l,
        const ctrs::array<typename state_t::value_type,3>& normal_r)
    {
        { t.calc_flux(ql, qr, normal_l, normal_r) } -> cvdf::ctrs::basic_array;
    };
    
    template <fluid_state::state_dependent_gas gas_t> struct totani_lr
    : public flux_base::flux_base_t<
        totani_lr,
        fluid_state::flux_t<dtype>,
        flux_input::flux_input_t<
            flux_input::left_right<cell_state<fluid_state::prim_t<dtype>>, cell_normal<fluid_state::prim_t<dtype>>>
            >
    {
        totani_lr(const gas_t& gas_in) {gas = &gas_in;}
        template <typename dtype> fluid_state::flux_t<dtype>
        calc_flux(
            const fluid_state::prim_t<dtype>& ql,
            const fluid_state::prim_t<dtype>& qr,
            const ctrs::array<dtype,3>& normal_l,
            const ctrs::array<dtype,3>& normal_r) const override final
        {
            fluid_state::flux_t<dtype> output;
            dtype un_l = normal_l[0]*ql.u()+normal_l[1]*ql.v()+normal_l[2]*ql.w();
            dtype un_r = normal_r[0]*qr.u()+normal_r[1]*qr.v()+normal_r[2]*qr.w();
            
            dtype rho_l = ql.p()/(gas->get_R(ql)*ql.T());
            dtype rho_r = qr.p()/(gas->get_R(qr)*qr.T());
            
            dtype e_l = ql.p()/(rho_l*(gas->get_gamma(ql)-1.0));
            dtype e_r = qr.p()/(rho_r*(gas->get_gamma(qr)-1.0));
            
            dtype c = 0.25*(rho_l+rho_r)*(un_l+un_r);
            output.continuity() = c;
            output.energy()     = 0.5*c*(e_l + e_r + ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w()) + 0.5*(un_l*qr.p() + un_r*ql.p());
            output.x_momentum() = 0.5*c*(ql.u()+qr.u()) + 0.5*(normal_l[0]*ql.p()+normal_r[0]*qr.p());
            output.y_momentum() = 0.5*c*(ql.v()+qr.v()) + 0.5*(normal_l[1]*ql.p()+normal_r[1]*qr.p());
            output.z_momentum() = 0.5*c*(ql.w()+qr.w()) + 0.5*(normal_l[2]*ql.p()+normal_r[2]*qr.p());
            return output;
        }
        
        const gas_t* gas;
    };
}