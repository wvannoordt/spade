#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/viscous_laws.h"
#include "navier-stokes/flux_input.h"

namespace cvdf::viscous
{
    template <viscous_laws::viscous_law vlaw_t> struct visc_lr
    {
        typedef typename vlaw_t::value_type dtype;
        typedef flux_input::flux_input_t
        <
            flux_input::left_right
            <
                flux_input::cell_info<>
            >,
            flux_input::face_info
            <
                flux_input::face_state<fluid_state::prim_t<dtype>>,
                flux_input::face_state_grad<ctrs::array<fluid_state::prim_t<dtype>, 3>>,
                flux_input::face_normal<ctrs::array<dtype, 3>>
            >
        > input_type;
        
        visc_lr(const vlaw_t& vlaw_in) { vlaw = &vlaw_in; }
        
        fluid_state::flux_t<dtype>
        calc_flux(const input_type& input) const
        {
            fluid_state::flux_t<dtype> output;
            // const auto& ql       = std::get<0>(input.cell_data.left.elements).data;
            // const auto& qr       = std::get<0>(input.cell_data.right.elements).data;
            // const auto& normal_l = std::get<1>(input.cell_data.left.elements).data;
            // const auto& normal_r = std::get<1>(input.cell_data.right.elements).data;
            
            output.continuity() = 0.0;
            output.energy()     = 0.0;
            output.x_momentum() = 0.0;
            output.y_momentum() = 0.0;
            output.z_momentum() = 0.0;
            return output;
        }
        
        const vlaw_t* vlaw;
    };
}