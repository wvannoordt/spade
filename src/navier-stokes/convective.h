#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/flux_base.h"

namespace cvdf::convective
{    
    template <fluid_state::state_dependent_gas gas_t> struct totani_lr
    : public flux_base::flux_base_t
    <
        totani_lr<gas_t>,
        fluid_state::flux_t<typename gas_t::value_type>,
        flux_base::flux_input_t
        <
            flux_base::left_right<>,
            flux_base::face_info<>
        >
    >
    {
        totani_lr(const gas_t& gas_in) {gas = &gas_in;}
        template <typename dtype> fluid_state::flux_t<dtype>
        calc_flux(const auto& input) const
        {
            fluid_state::flux_t<dtype> output;
            // const auto& ql = input.cell_data.left.get<0>.data;
            // const auto& qr = input.cell_data.right.get<0>.data;
            // const auto& normal_l = input.cell_data.left.get<1>.data;
            // const auto& normal_r = input.cell_data.right.get<1>.data;
            
            // dtype un_l = normal_l[0]*ql.u()+normal_l[1]*ql.v()+normal_l[2]*ql.w();
            // dtype un_r = normal_r[0]*qr.u()+normal_r[1]*qr.v()+normal_r[2]*qr.w();
            // 
            // dtype rho_l = ql.p()/(gas->get_R(ql)*ql.T());
            // dtype rho_r = qr.p()/(gas->get_R(qr)*qr.T());
            // 
            // dtype e_l = ql.p()/(rho_l*(gas->get_gamma(ql)-1.0));
            // dtype e_r = qr.p()/(rho_r*(gas->get_gamma(qr)-1.0));
            // 
            // dtype c = 0.25*(rho_l+rho_r)*(un_l+un_r);
            // output.continuity() = c;
            // output.energy()     = 0.5*c*(e_l + e_r + ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w()) + 0.5*(un_l*qr.p() + un_r*ql.p());
            // output.x_momentum() = 0.5*c*(ql.u()+qr.u()) + 0.5*(normal_l[0]*ql.p()+normal_r[0]*qr.p());
            // output.y_momentum() = 0.5*c*(ql.v()+qr.v()) + 0.5*(normal_l[1]*ql.p()+normal_r[1]*qr.p());
            // output.z_momentum() = 0.5*c*(ql.w()+qr.w()) + 0.5*(normal_l[2]*ql.p()+normal_r[2]*qr.p());
            return output;
        }
        
        const gas_t* gas;
    };
}