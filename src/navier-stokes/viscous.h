#pragma once

#include "core/config.h"
#include "core/ctrs.h"

#include "omni/omni.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/viscous_laws.h"

namespace spade::viscous
{
    template <typename vlaw_t, typename gas_model_t> struct visc_lr
    {
        using float_t       = vlaw_t::value_type;
        using output_type   = fluid_state::flux_t<float_t>;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::gradient, omni::info::metric>;
        using sub_info_type = omni::info_union<typename vlaw_t::info_type, typename gas_model_t::info_type>;
        using info_type     = omni::info_union<own_info_type, sub_info_type>;
        using omni_type     = omni::prefab::face_mono_t<info_type>;
        
        visc_lr(const vlaw_t& vlaw_in, const gas_model_t& gas_in) : vlaw{vlaw_in}, gas{gas_in} {}
        
        _sp_hybrid output_type operator() (const auto& input) const
        {
            // don't forget the negative signs
            output_type output;
            const auto& q_face      = omni::access<omni::info::value>   (input.face(0_c));
            const auto& face_grad   = omni::access<omni::info::gradient>(input.face(0_c));
            const auto& n           = omni::access<omni::info::metric>  (input.face(0_c));
            
            const auto& visc        = vlaw.get_visc(input.face(0_c));
            const auto& visc2       = vlaw.get_beta(input.face(0_c));
            const auto& div         = face_grad[0].u() + face_grad[1].v() + face_grad[2].w();
            linear_algebra::dense_mat<float_t, 3> tau;
            ctrs::array<float_t, 3> ht;
            
            tau(0,0) = visc*(face_grad[0].u() + face_grad[0].u()) + visc2*div; //tau_xx
            tau(0,1) = visc*(face_grad[0].v() + face_grad[1].u()); //tau_xy
            tau(0,2) = visc*(face_grad[0].w() + face_grad[2].u()); //tau_xz
            tau(1,0) = visc*(face_grad[1].u() + face_grad[0].v()); //tau_yx
            tau(1,1) = visc*(face_grad[1].v() + face_grad[1].v()) + visc2*div; //tau_yy
            tau(1,2) = visc*(face_grad[1].w() + face_grad[2].v()); //tau_yz
            tau(2,0) = visc*(face_grad[2].u() + face_grad[0].w()); //tau_zx
            tau(2,1) = visc*(face_grad[2].v() + face_grad[1].w()); //tau_zy
            tau(2,2) = visc*(face_grad[2].w() + face_grad[2].w()) + visc2*div; //tau_zz
  
            const auto gam       = gas.get_gamma(input.face(0_c));
            const auto rgas      = gas.get_R(input.face(0_c));
            const auto spec_heat = gam*rgas/(gam-float_t(1.0));
            const auto dfs       = vlaw.get_diffuse(input.face(0_c));
            const auto cond      = spec_heat*dfs;


            ht[0] = q_face.u()*tau(0,0)+q_face.v()*tau(0,1)+q_face.w()*tau(0,2) + cond*face_grad[0].T();
            ht[1] = q_face.u()*tau(1,0)+q_face.v()*tau(1,1)+q_face.w()*tau(1,2) + cond*face_grad[1].T();
            ht[2] = q_face.u()*tau(2,0)+q_face.v()*tau(2,1)+q_face.w()*tau(2,2) + cond*face_grad[2].T();
            output.continuity() = 0.0;
            output.energy()     = -(n[0]*ht[0]+n[1]*ht[1]+n[2]*ht[2]);
            output.x_momentum() = -(n[0]*tau(0,0)+n[1]*tau(0,1)+n[2]*tau(0,2));
            output.y_momentum() = -(n[0]*tau(1,0)+n[1]*tau(1,1)+n[2]*tau(1,2));
            output.z_momentum() = -(n[0]*tau(2,0)+n[1]*tau(2,1)+n[2]*tau(2,2));
            return output;
        }
        
        const vlaw_t vlaw;
        const gas_model_t gas;
    };
}