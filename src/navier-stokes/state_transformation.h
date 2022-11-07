#pragma once

#include "core/algs.h"

#include "navier-stokes/fluid_state.h"

namespace spade::fluid_state
{
    namespace detail
    {
        template <typename gas_t, fluid_state::is_state_type from_t, fluid_state::is_state_type to_t>
        requires (fluid_state::state_convertible<from_t, to_t, gas_t>)
        struct mono_state_converstion_t
        {
            const gas_t* gas;
            typedef from_t arg_type;
            mono_state_converstion_t(const gas_t& gas_in) {gas = &gas_in;}
            to_t operator () (const from_t& q) const
            {
                to_t w;
                spade::fluid_state::convert_state(q, w, *gas);
                return w;
            };
        };
    }
    
    template <typename gas_t, fluid_state::is_state_type forward_t>
    struct state_transform_t
    {
        const gas_t* gas;
        
        //prim = inverse
        //cons = forward
        
        state_transform_t(const forward_t& to_state, const gas_t& gas_in) { gas = &gas_in; }
        
        template <typename array_t>
        requires (fluid_state::is_state_type<typename array_t::alias_type> && fluid_state::state_convertible<forward_t, typename array_t::alias_type, gas_t>)
        void transform_forward (array_t& q) const
        {
            using inverse_t = typename array_t::alias_type;
            using i2f_t = detail::mono_state_converstion_t<gas_t, inverse_t, forward_t>;
            spade::algs::transform_inplace(q, i2f_t(*gas));
        }
        
        template <typename array_t>
        requires (fluid_state::is_state_type<typename array_t::alias_type> && fluid_state::state_convertible<forward_t, typename array_t::alias_type, gas_t>)
        void transform_inverse (array_t& q) const
        {
            using inverse_t = typename array_t::alias_type;
            using f2i_t = detail::mono_state_converstion_t<gas_t, forward_t, inverse_t>;
            spade::algs::transform_inplace(q, f2i_t(*gas));
        }
    };
}