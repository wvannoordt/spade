#pragma once

#include "core/algs.h"
#include "core/unsafe_cast.h"
#include "omni/omni.h"

#include "navier-stokes/fluid_state.h"

namespace spade::fluid_state
{
    namespace detail
    {
        template <typename array_t, typename gas_t, fluid_state::is_state_type from_t, fluid_state::is_state_type to_t>
        requires (fluid_state::state_convertible<from_t, to_t, gas_t>)
        struct mono_state_converstion_t
        {
            constexpr static grid::array_centering centr = array_t::centering_type();
            using native_t  = array_t::alias_type;
            using omni_type = omni::prefab::mono_t<centr, omni::info_list_t<omni::info::value>>;
            const gas_t gas;
            mono_state_converstion_t(const gas_t& gas_in) : gas{gas_in} {}
            _sp_hybrid native_t operator () (const auto& input) const
            {
                // This is the worst thing in the whole universe,
                // if you're an employer and you are seeing this then please realize that
                // I was young and hopelessly naive when I wrote this
                auto state = omni::access<omni::info::value>(input.template seek_element<centr>(0_c));
                from_t& q = algs::unsafe_cast<from_t>(state);
                to_t w;
                spade::fluid_state::convert_state(q, w, gas);
                return algs::unsafe_cast<native_t>(w);
            };
        };
    }
    
    template <typename gas_t, fluid_state::is_state_type forward_t>
    struct state_transform_t
    {
        const gas_t& gas;
        const grid::exchange_inclusion_e g_config;
        
        //prim = inverse
        //cons = forward
        
        state_transform_t(const forward_t&, const gas_t& gas_in)                                         : gas{gas_in}, g_config{grid::include_exchanges} {}
        state_transform_t(const forward_t&, const gas_t& gas_in, const grid::exchange_inclusion_e& g_in) : gas{gas_in}, g_config{g_in} {}
        
        template <typename array_t>
        requires (fluid_state::is_state_type<typename array_t::alias_type> && fluid_state::state_convertible<forward_t, typename array_t::alias_type, gas_t>)
        void transform_forward (array_t& q) const
        {
            using inverse_t = typename array_t::alias_type;
            using i2f_t = detail::mono_state_converstion_t<array_t, gas_t, inverse_t, forward_t>;
            spade::algs::transform_inplace(q, i2f_t(gas), g_config);
        }
        
        template <typename array_t>
        requires (fluid_state::is_state_type<typename array_t::alias_type> && fluid_state::state_convertible<forward_t, typename array_t::alias_type, gas_t>)
        void transform_inverse (array_t& q) const
        {
            using inverse_t = typename array_t::alias_type;
            using f2i_t = detail::mono_state_converstion_t<array_t, gas_t, forward_t, inverse_t>;
            spade::algs::transform_inplace(q, f2i_t(gas), g_config);
        }
    };
}