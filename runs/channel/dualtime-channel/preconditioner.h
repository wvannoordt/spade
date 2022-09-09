#pragma once
#include "spade.h"
template <typename gas_t, typename array_t, typename beta_t>
struct preconditioner_t
{
    const gas_t* gas;
    beta_t beta;
    preconditioner_t(const gas_t& gas_in, const array_t& q, const beta_t& beta_in)
    {
        gas  = &gas_in;
        beta = beta_in;
    }
    
    void apply_scale_factor(array_t& q, const beta_t& factor) const
    {
        struct inner_trans_t
        {
            beta_t fac;
            typedef cons_t arg_type;
            cons_t operator() (const arg_type& c)
            {
                prim_t p;
                cons_t w;
                spade::fluid_state::convert_state(c, p, *gas);
                p.p() *= fac;
                spade::fluid_state::convert_state(p, w, *gas);
                return w;
            }
        } inner_trans{.fac=factor};
        spade::algs::transform_inplace(q, inner_trans);
    }
    
    void transform_forward(array_t& q) const
    {
        apply_scale_factor(q, beta);
    }
    
    void transform_inverse(array_t& q) const
    {
        apply_scale_factor(q, 1.0/beta);
    }
};