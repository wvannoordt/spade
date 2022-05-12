#pragma once
#include <concepts>

#include "core/ctrs.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/viscous_laws.h"

namespace cvdf::viscous
{
    template <viscous_laws::viscous_law vlaw_t> struct visc_lr
    {
        visc_lr(const vlaw_t& vlaw_in)
        {
            vlaw = &vlaw_in;
        }
        
        const vlaw_t* vlaw;
    };
}