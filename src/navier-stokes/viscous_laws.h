#pragma once

#include "navier-stokes/fluid_state.h"
#include <concepts>


namespace cvdf::viscous_laws
{
    template <typename dtype> struct constant_viscosity_t
    {
        constant_viscosity_t(dtype visc_in) { this->visc = visc_in; this->beta = 0.66666666667*this->visc; }
        constant_viscosity_t(void) { this->visc = dtype(); this->beta = 0.66666666667*this->visc; }
        
        // template <prim_state_type state_type>dtype calc_visc(const state_type& q)
        template <class state_type>dtype calc_visc(const state_type& q)
        {
            return visc;
        }
        
        dtype visc;
        dtype prandtl;
        dtype beta;
    };
}