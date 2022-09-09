#pragma once

#include "spade.h"
#include "typedef.h"

template <typename array_t> struct trans_t
{
    using gas_t = spade::fluid_state::perfect_gas_t<real_t>;
    const gas_t* gas;
    
    struct p2c_t
    {
        const gas_t* gas;
        typedef prim_t arg_type;
        p2c_t(const gas_t& gas_in) {gas = &gas_in;}
        cons_t operator () (const prim_t& q) const
        {
            cons_t w;
            spade::fluid_state::convert_state(q, w, *gas);
            return w;
        };
    };
    struct c2p_t
    {
        const gas_t* gas;
        typedef cons_t arg_type;
        c2p_t(const gas_t& gas_in) {gas = &gas_in;}
        prim_t operator () (const cons_t& w) const
        {
            prim_t q;
            spade::fluid_state::convert_state(w, q, *gas);
            return q;
        }
    };
    
    trans_t(const gas_t& gas_in, const array_t& q) { gas = &gas_in; }
    
    void transform_forward (array_t& q) const { spade::algs::transform_inplace(q, p2c_t(*gas)); }
    void transform_inverse (array_t& q) const { spade::algs::transform_inplace(q, c2p_t(*gas)); }
};