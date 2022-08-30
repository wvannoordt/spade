#pragma once

#include <concepts>



namespace spade::time_integration
{
    template <typename T> concept time_int_scheme = requires(T t)
    {
        T::num_substeps();
        t.advance();
    };
    
    struct identity_transform_t { };
    
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
        
        static std::size_t num_substeps(void) {return 2;}
        
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
    
    // template <
    //     typename var_state_t,
    //     typename rhs_state_t,
    //     typename time_state_t,
    //     typename rhs_calc_t,
    //     typename rk_table_t=rk_tables::rk4_t,
    //     typename state_trans_t=identity_transform_t,
    //     typename inv_state_trans_t=identity_transform_t
    //     >
    // struct explicit_rk
    // {
    //     explicit_rk(
    //         const var_state_t& q_in,
    //         const rhs_state_t& rhs_in,
    //         const time_state_t& t0_in,
    //         const time_state_t& dt_in,
    //         const rhs_calc_t& rhs_calc_in,
    //         // const rk_table_t& rk_table_in = rk_table_t(),
    //         const state_trans_t& state_trans_in = identity_transform_t(),
    //         const inv_state_trans_t& inv_state_trans_in = identity_transform_t())
    //     {
    // 
    //     }
    // 
    //     // rk_table_t table;
    // };
    // 
    // template <
    //     typename var_state_t,
    //     typename rhs_state_t,
    //     typename time_state_t,
    //     typename rhs_calc_t,
    //     typename rk_table_t=rk_tables::rk4_t,
    //     typename state_trans_t=identity_transform_t,
    //     typename inv_state_trans_t=identity_transform_t
    //     >
    // https://en.wikipedia.org/wiki/Backward_differentiation_formula
    // struct implicit_bdf
    // {
    //     explicit_rk(
    //         const var_state_t& q_in,
    //         const rhs_state_t& rhs_in,
    //         const time_state_t& t0_in,
    //         const time_state_t& dt_in,
    //         const rhs_calc_t& rhs_calc_in,
    //         // const rk_table_t& rk_table_in = rk_table_t(),
    //         const state_trans_t& state_trans_in = identity_transform_t(),
    //         const inv_state_trans_t& inv_state_trans_in = identity_transform_t())
    //     {
    // 
    //     }
    // 
    //     // rk_table_t table;
    // };
    
    template
    <
        const std::size_t order=2,
        typename var_state_t,
        typename residual_state_t,
        typename time_state_t,
        typename rhs_calc_t,
        typename implicit_solve_t,
        typename state_trans_t=identity_transform_t,
        typename inv_state_trans_t=identity_transform_t
        
    > struct bdf_t
    {
        
    };
}