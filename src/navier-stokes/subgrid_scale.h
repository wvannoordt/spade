#pragma once

#include "omni/omni.h"

namespace spade::subgrid_scale
{
    template <typename float_t> struct wale_t
    {
        float_t cw, delta;
        wale_t(const float_t& cw_in, const float_t& delta_in) : cw{cw_in}, delta{delta_in} {}
        using omni_type = omni::prefab::mono_t<grid::agno_centered, omni::info::gradient>;
        float_t get_mut (const auto& input) const
        {
            const auto& grad = omni::access<omni::info::gradient>(input.root());
            return grad[0].u() + grad[1].v() + grad[2].w();
        }
    };
}