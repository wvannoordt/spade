#pragma once

#include "core/config.h"
#include "core/ctrs.h"
#include "fetch/fetch.h"

#include "navier-stokes/fluid_state.h"

namespace spade::convective
{
    template <typename dtype> using standard_lr_input_type = fetch::face_fetch_t
    <
        fetch::left_right
        <
            fetch::cell_info
            <
                fetch::cell_state<fluid_state::prim_t<dtype>>,
                fetch::cell_normal<ctrs::array<dtype, 3>>
            >
        >,
        fetch::face_info<>
    >;
        
    template <fluid_state::state_dependent_gas gas_t> struct totani_lr
    {
        typedef typename gas_t::value_type dtype;
        typedef fluid_state::flux_t<dtype> output_type;
        typedef standard_lr_input_type<dtype> input_type;
        
        totani_lr(const gas_t& gas_in) {gas = &gas_in;}
        
        
        output_type calc_flux(const input_type& input) const
        {
            output_type output;
            const auto& ql       = std::get<0>(input.cell_data.left.elements).data;
            const auto& qr       = std::get<0>(input.cell_data.right.elements).data;
            const auto& normal_l = std::get<1>(input.cell_data.left.elements).data;
            const auto& normal_r = std::get<1>(input.cell_data.right.elements).data;
            
            dtype un_l = normal_l[0]*ql.u()+normal_l[1]*ql.v()+normal_l[2]*ql.w();
            dtype un_r = normal_r[0]*qr.u()+normal_r[1]*qr.v()+normal_r[2]*qr.w();
            
            dtype rho_l = ql.p()/(gas->get_R(ql)*ql.T());
            dtype rho_r = qr.p()/(gas->get_R(qr)*qr.T());
            
            dtype e_l = ql.p()/(rho_l*(gas->get_gamma(ql)-dtype(1.0)));
            dtype e_r = qr.p()/(rho_r*(gas->get_gamma(qr)-dtype(1.0)));
            
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
    
    template <fluid_state::state_dependent_gas gas_t> struct pressure_diss_lr
    {
        typedef typename gas_t::value_type dtype;
        typedef fluid_state::flux_t<dtype> output_type;
        typedef standard_lr_input_type<dtype> input_type;
        
        pressure_diss_lr(const gas_t& gas_in, const dtype& eps_p_in) {gas = &gas_in; eps_p = eps_p_in; eps_T = eps_p_in;}
        pressure_diss_lr(const gas_t& gas_in, const dtype& eps_p_in, const dtype& eps_T_in) {gas = &gas_in; eps_p = eps_p_in; eps_T = eps_T_in;}
        
        output_type calc_flux(const input_type& input) const
        {
            output_type output;
            const auto& ql       = std::get<0>(input.cell_data.left.elements).data;
            const auto& qr       = std::get<0>(input.cell_data.right.elements).data;
            const auto& normal_l = std::get<1>(input.cell_data.left.elements).data;
            const auto& normal_r = std::get<1>(input.cell_data.right.elements).data;
            
            const dtype delta_p = ql.p()-qr.p();
            const dtype delta_T = ql.T()-qr.T();
            const dtype pavg = 0.5*(ql.p() + qr.p());
            const dtype Tavg = 0.5*(ql.T() + qr.T());
            const dtype inv_RT = 1.0/(0.5*(gas->get_R(ql)+gas->get_R(qr))*Tavg);
            const dtype p_inv_rt2 = pavg*inv_RT/Tavg;
            const dtype gamma =  0.5*(gas->get_gamma(ql)+gas->get_gamma(qr));
            const dtype k = 0.5*(ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w());
            const dtype u = 0.5*(ql.u()+qr.u());
            const dtype v = 0.5*(ql.v()+qr.v());
            const dtype w = 0.5*(ql.w()+qr.w());
            output.continuity() = eps_p*delta_p*inv_RT;
            output.energy()     = eps_p*delta_p*(k*inv_RT + 1.0/gamma);
            output.x_momentum() = eps_p*delta_p*u*inv_RT;
            output.y_momentum() = eps_p*delta_p*v*inv_RT;
            output.z_momentum() = eps_p*delta_p*w*inv_RT;
            
            output.continuity() -= eps_T*delta_T*p_inv_rt2;
            output.energy()     -= eps_T*delta_T*p_inv_rt2*(0.5*k);
            output.x_momentum() -= eps_T*delta_T*p_inv_rt2*u;
            output.y_momentum() -= eps_T*delta_T*p_inv_rt2*v;
            output.z_momentum() -= eps_T*delta_T*p_inv_rt2*w;
            
            return output;
        }
        
        const gas_t* gas;
        dtype eps_p, eps_T;
    };
    
    template <fluid_state::state_dependent_gas gas_t> struct weno_3
    {
        typedef typename gas_t::value_type dtype;
        typedef fluid_state::flux_t<dtype> output_type;
        typedef fetch::face_fetch_t
        <
            fetch::flux_line
            <
                4,
                fetch::cell_info
                <
                    fetch::cell_state<fluid_state::prim_t<dtype>>,
                    fetch::cell_normal<ctrs::array<dtype, 3>>
                >
            >,
            fetch::face_info<>
        > input_type;
        
        weno_3(const gas_t& gas_in) { gas = &gas_in; }
        
        output_type calc_flux(const input_type& input) const
        {
            fluid_state::flux_t<dtype> output;
            const auto& q0       = std::get<0>(input.cell_data.stencil[0].elements).data;
            const auto& q1       = std::get<0>(input.cell_data.stencil[1].elements).data;
            const auto& q2       = std::get<0>(input.cell_data.stencil[2].elements).data;
            const auto& q3       = std::get<0>(input.cell_data.stencil[3].elements).data;
            const auto& n0       = std::get<1>(input.cell_data.stencil[0].elements).data;
            const auto& n1       = std::get<1>(input.cell_data.stencil[1].elements).data;
            const auto& n2       = std::get<1>(input.cell_data.stencil[2].elements).data;
            const auto& n3       = std::get<1>(input.cell_data.stencil[3].elements).data;
            
            fluid_state::flux_t<dtype> f0_u, f0_d, f1_u, f1_d, f2_u, f2_d, f3_u, f3_d;
            const auto flux_q = [&](
                const fluid_state::prim_t<dtype>& q,
                const ctrs::array<dtype, 3>& nv,
                fluid_state::flux_t<dtype>& f_u,
                fluid_state::flux_t<dtype>& f_d) -> void
            {
                const dtype u_n  = nv[0]*q.u() + nv[1]*q.v() + nv[2]*q.w();
                const dtype rho  = q.p()/(gas->get_R(q)*q.T());
                const dtype enth = 0.5*(q.u()*q.u()+q.v()*q.v()+q.w()*q.w()) + q.p()/rho + q.p()/(rho*(gas->get_gamma(q)-1.0));
                f_u.continuity() = 0.5*rho*u_n;
                f_u.energy()     = 0.5*rho*u_n*enth;
                f_u.x_momentum() = 0.5*rho*q.u()*u_n + 0.5*q.p()*nv[0];
                f_u.y_momentum() = 0.5*rho*q.v()*u_n + 0.5*q.p()*nv[1];
                f_u.z_momentum() = 0.5*rho*q.w()*u_n + 0.5*q.p()*nv[2];
                f_d = f_u;
                const real_t sigma = abs(u_n) + sqrt(gas->get_gamma(q)*gas->get_R(q)*q.T());
                f_d.continuity() -= sigma*rho;
                f_d.energy()     -= sigma*enth;
                f_d.x_momentum() -= sigma*rho*q.u();
                f_d.y_momentum() -= sigma*rho*q.v();
                f_d.z_momentum() -= sigma*rho*q.w();
                f_u.continuity() += sigma*rho;
                f_u.energy()     += sigma*enth;
                f_u.x_momentum() += sigma*rho*q.u();
                f_u.y_momentum() += sigma*rho*q.v();
                f_u.z_momentum() += sigma*rho*q.w();
            };
            flux_q(q0, n0, f0_u, f0_d);
            flux_q(q1, n1, f1_u, f1_d);
            flux_q(q2, n2, f2_u, f2_d);
            flux_q(q3, n3, f3_u, f3_d);
            
            for (std::size_t k = 0; k < output.size(); ++k)
            {
                //upwind
                const dtype r0  = -0.5*f0_u[k]+1.5*f1_u[k]; //sten0
                const dtype r1  =  0.5*f1_u[k]+0.5*f2_u[k]; //sten1
                
                //downwind
                const dtype r2  =  0.5*f1_d[k]+0.5*f2_d[k]; //sten2
                const dtype r3  =  1.5*f2_d[k]-0.5*f3_d[k]; //sten3
                
                const dtype be0 = (f0_u[k]-f1_u[k])*(f0_u[k]-f1_u[k]);
                const dtype be1 = (f1_u[k]-f2_u[k])*(f1_u[k]-f2_u[k]);
                
                const dtype be2 = (f1_d[k]-f2_d[k])*(f1_d[k]-f2_d[k]);
                const dtype be3 = (f2_d[k]-f3_d[k])*(f2_d[k]-f3_d[k]);
                
                const dtype eps = 1e-9;
                const dtype a0  = (1.0/3.0)/((be0+eps)*(be0+eps));
                const dtype a1  = (2.0/3.0)/((be1+eps)*(be1+eps));
                const dtype a2  = (2.0/3.0)/((be2+eps)*(be2+eps));
                const dtype a3  = (1.0/3.0)/((be3+eps)*(be3+eps));
                
                const dtype w0 = a0/(a0+a1);
                const dtype w1 = a1/(a0+a1);
                const dtype w2 = a2/(a2+a3);
                const dtype w3 = a3/(a2+a3);
                
                output[k] = w0*r0 + w1*r1 + w2*r2 + w3*r3;
                // output[k] = (1.0/3.0)*r0 + (2.0/3.0)*r1 + (2.0/3.0)*r2 + (1.0/3.0)*r3; //(mms)
            }
            
            return output;
        }
        
        const gas_t* gas;
    };
}
