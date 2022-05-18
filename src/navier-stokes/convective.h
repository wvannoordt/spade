#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/flux_input.h"

namespace cvdf::convective
{    
    template <fluid_state::state_dependent_gas gas_t> struct totani_lr
    {
        typedef flux_input::flux_input_t
        <
            flux_input::left_right
            <
                flux_input::cell_info
                <
                    flux_input::cell_state<fluid_state::prim_t<typename gas_t::value_type>>,
                    flux_input::cell_normal<ctrs::array<typename gas_t::value_type, 3>>
                >
            >,
            flux_input::face_info<>
        > input_type;
        
        typedef typename gas_t::value_type dtype;
        
        totani_lr(const gas_t& gas_in) {gas = &gas_in;}
        
        fluid_state::flux_t<dtype>
        calc_flux(const input_type& input) const
        {
            fluid_state::flux_t<dtype> output;
            const auto& ql       = std::get<0>(input.cell_data.left.elements).data;
            const auto& qr       = std::get<0>(input.cell_data.right.elements).data;
            const auto& normal_l = std::get<1>(input.cell_data.left.elements).data;
            const auto& normal_r = std::get<1>(input.cell_data.right.elements).data;
            
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