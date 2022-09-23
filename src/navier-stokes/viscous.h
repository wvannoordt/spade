#pragma once
#include <concepts>

#include "core/ctrs.h"
#include "core/fetch.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/viscous_laws.h"

namespace spade::viscous
{
    template <viscous_laws::viscous_law vlaw_t> struct visc_lr
    {
        typedef typename vlaw_t::value_type dtype;
        typedef fluid_state::flux_t<dtype> output_type;
        typedef fetch::face_fetch_t
        <
            fetch::left_right
            <
                fetch::cell_info<>
            >,
            fetch::face_info
            <
                fetch::face_state<fluid_state::prim_t<dtype>>,
                fetch::face_state_grad<ctrs::array<fluid_state::prim_t<dtype>, 3>>,
                fetch::face_normal<ctrs::array<dtype, 3>>
            >
        > input_type;
        
        visc_lr(const vlaw_t& vlaw_in) { vlaw = &vlaw_in; }
        
        output_type calc_flux(const input_type& input) const
        {
            // don't forget the negative signs
            fluid_state::flux_t<dtype> output;
            const auto& q_face      = std::get<0>(input.face_data.elements).data;
            const auto& face_grad   = std::get<1>(input.face_data.elements).data;
            const auto& n           = std::get<2>(input.face_data.elements).data;
            const auto& visc        = vlaw->get_visc(q_face);
            const auto& visc2       = vlaw->get_beta(q_face);
            const auto& div         = face_grad[0].u() + face_grad[1].v() + face_grad[2].w();
            linear_algebra::dense_mat<dtype, 3> tau;
            ctrs::array<dtype, 3> ht;
            
            tau(0,0) = visc*(face_grad[0].u() + face_grad[0].u()) + visc2*div; //tau_xx
            tau(0,1) = visc*(face_grad[0].v() + face_grad[1].u()); //tau_xy
            tau(0,2) = visc*(face_grad[0].w() + face_grad[2].u()); //tau_xz
            tau(1,0) = visc*(face_grad[1].u() + face_grad[0].v()); //tau_yx
            tau(1,1) = visc*(face_grad[1].v() + face_grad[1].v()) + visc2*div; //tau_yy
            tau(1,2) = visc*(face_grad[1].w() + face_grad[2].v()); //tau_yz
            tau(2,0) = visc*(face_grad[2].u() + face_grad[0].w()); //tau_zx
            tau(2,1) = visc*(face_grad[2].v() + face_grad[1].w()); //tau_zy
            tau(2,2) = visc*(face_grad[2].w() + face_grad[2].w()) + visc2*div; //tau_zz
            
            ht[0] = q_face.u()*tau(0,0)+q_face.v()*tau(0,1)+q_face.w()*tau(0,2) - vlaw->get_conductivity(q_face)*face_grad[0].T();
            ht[1] = q_face.u()*tau(1,0)+q_face.v()*tau(1,1)+q_face.w()*tau(1,2) - vlaw->get_conductivity(q_face)*face_grad[1].T();
            ht[2] = q_face.u()*tau(2,0)+q_face.v()*tau(2,1)+q_face.w()*tau(2,2) - vlaw->get_conductivity(q_face)*face_grad[2].T();
            output.continuity() = 0.0;
            output.energy()     = -(n[0]*ht[0]+n[1]*ht[1]+n[2]*ht[2]);
            output.x_momentum() = -(n[0]*tau(0,0)+n[1]*tau(0,1)+n[2]*tau(0,2));
            output.y_momentum() = -(n[0]*tau(1,0)+n[1]*tau(1,1)+n[2]*tau(1,2));
            output.z_momentum() = -(n[0]*tau(2,0)+n[1]*tau(2,1)+n[2]*tau(2,2));
            return output;
        }
        
        const vlaw_t* vlaw;
    };
}