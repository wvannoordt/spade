#pragma once

namespace cvdf::time_integration
{
    template <typename state_trans_t> struct identity_transform_t
    {
        void operator() (state_trans_t& state) const {}
    };
    template <
        typename var_state_t,
        typename rhs_state_t,
        typename time_state_t,
        typename rhs_calc_t,
        typename state_trans_t,
        typename inv_state_trans_t
        >
    struct rk2
    {
        time_state_t t;
        time_state_t dt;
        var_state_t* q;
        rhs_state_t* k0;
        rhs_state_t* k1;
        const rhs_calc_t*  rhs_oper;
        const state_trans_t* state_trans;
        const inv_state_trans_t* inv_state_trans;
        rk2(
            var_state_t& q_in,
            rhs_state_t& k0_in,
            rhs_state_t& k1_in,
            const time_state_t& t0_in,
            const time_state_t& dt_in,
            const rhs_calc_t& rhs_oper_in,
            const state_trans_t& trans_in,
            const inv_state_trans_t& inv_trans_in)
        {
            dt = dt_in;
            q  = &q_in;
            t  = t0_in;
            k0 = &k0_in;
            k1 = &k1_in;
            rhs_oper = &rhs_oper_in;
            state_trans = &trans_in;
            inv_state_trans = &inv_trans_in;
        }
        
        void advance()
        {
            (*rhs_oper)((*k0), (*q), t);
            (*k0) *= (0.5*dt);
            (*state_trans)((*q));
            (*q) += (*k0);
            (*inv_state_trans)((*q));
            t += 0.5*dt;
            (*rhs_oper)((*k1), (*q), t);
            (*k1) *= (dt);
            (*state_trans)((*q));
            (*q) -= (*k0);
            (*q) += (*k1);
            (*inv_state_trans)((*q));
            t += 0.5*dt;
        }
        
        time_state_t& time(void) {return  t;}
        var_state_t&  soln(void) {return *q;}
    };
}