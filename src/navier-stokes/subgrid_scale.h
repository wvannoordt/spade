#pragma once

#include "omni/omni.h"

namespace spade::subgrid_scale
{
    template <typename derived_t, omni::info::is_omni_info... infos_t>
    struct subgrid_interface_t
    {
        derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = omni::info_list_t<infos_t...>;

        auto get_mut(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_mut(args...);}, input);
        }
    };

    template <typename float_t> struct wale_t
    : public subgrid_interface_t<wale_t<float_t>, omni::info::gradient>
    {
        using base_t = subgrid_interface_t<wale_t<float_t>, omni::info::gradient>;
        using base_t::get_mut;
        using base_t::info_type;

        float_t cw, delta;
        wale_t(const float_t& cw_in, const float_t& delta_in) : cw{cw_in}, delta{delta_in} {}

        float_t get_mut (const ctrs::array<fluid_state::prim_t<float_t>, 3>& grad) const
        {
            return grad[0].u() + grad[1].v() + grad[2].w();
        }
    };
}