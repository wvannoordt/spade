#pragma once

#include <concepts>
#include "core/finite_diff.h"
#include "core/ctrs.h"

namespace spade::time_integration
{
    enum var_cascade_policy
    {
        var_cascade_assign,
        var_cascade_inplace
    };
    
    namespace detail
    {
        template <const std::size_t ar_size, const var_cascade_policy policy> struct cascade_offset_t {};
        template <const std::size_t ar_size> struct cascade_offset_t<ar_size, var_cascade_assign> {};
        template <const std::size_t ar_size> struct cascade_offset_t<ar_size, var_cascade_inplace>
        {
            std::size_t value = 0;
            cascade_offset_t& operator ++()
            {
                value++;
                value %= ar_size;
                return *this;
            }
        };
        
        template <typename cas_t>
        requires (cas_t::get_policy()==var_cascade_inplace)
        static std::size_t get_index(const cas_t& cas, const std::size_t& idx)
        {
            return (idx+cas.offset.value) % cas.size();
        }
        
        template <typename cas_t>
        requires (cas_t::get_policy()==var_cascade_assign)
        static std::size_t get_index(const cas_t& cas, const std::size_t& idx)
        {
            return idx;
        }
        
        template <typename cas_t, typename state_t>
        requires (cas_t::get_policy()==var_cascade_inplace)
        static void cascade(cas_t& cas, const state_t& q)
        {
            ++cas.offset;
            cas.instances[cas.index(cas.size()-1)] = q;
        }
        
        template <typename cas_t, typename state_t>
        requires (cas_t::get_policy()==var_cascade_assign)
        static void cascade(cas_t& cas, const state_t& q)
        {
            for (auto i: range(0,cas.size()-1))
            {
                cas.instances[cas.index(i)] = cas.instances[cas.index(i+1)];
            }
            cas.instances[cas.index(cas.size()-1)] = q;
        }
        
        template <typename state_t, const std::size_t ar_size, const var_cascade_policy policy> struct cascade_t
        {
            constexpr static var_cascade_policy get_policy() {return policy;}
            ctrs::array<state_t, ar_size> instances;
            cascade_offset_t<ar_size, policy> offset;
            constexpr static std::size_t size() {return ar_size;}
            std::size_t index(const std::size_t& idx) const { return get_index(*this, idx); }
            const state_t& operator [](const std::size_t& idx) const {return instances[this->index(idx)];}
            state_t& operator [](const std::size_t& idx) {return instances[this->index(idx)];}
            void push(const state_t& q)
            {
                cascade(*this, q);
            }
        };
        
        template <typename T> concept has_fundamental_type = requires(T t) { T::fundamental_type; };
        template <typename T> concept has_value_type = requires(T t) { T::value_type; };
        
        template <typename default_t, typename query_t> struct get_fundamental_type_if_avail
        {
            typedef default_t type;
        };
        
        template <typename default_t, has_fundamental_type query_t> struct get_fundamental_type_if_avail<default_t, query_t>
        {
            typedef typename query_t::fundamental_type type;
        };
        
        template <typename default_t, has_value_type query_t> struct get_fundamental_type_if_avail<default_t, query_t>
        {
            typedef typename query_t::value_type type;
        };
    }
    
    template <typename T> concept time_int_scheme = requires(T t)
    {
        T::num_substeps();
        t.advance();
    };
    
    struct identity_transform_t
    {
        template <typename input_t> void operator()(const input_t& input) const {}
    };
    
    static identity_transform_t identity_transform;
    
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
        
        time_state_t& time(void)     {return  t;}
        var_state_t&  solution(void) {return *q;}
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
    
    const static int default_bdf_order = 2;
    template
    <
        typename var_state_t,
        typename residual_state_t,
        typename time_state_t,
        typename rhs_calc_t,
        typename implicit_solve_t,
        const int order=default_bdf_order,
        typename state_trans_t=identity_transform_t,
        typename inv_state_trans_t=identity_transform_t,
        const var_cascade_policy policy=var_cascade_inplace
    > struct bdf_t
    {
        static_assert(order <= 6, "BDF non-zero-stable for order larger than 6!");
        typedef typename detail::get_fundamental_type_if_avail<time_state_t, var_state_t>::type coeff_t;
        
        var_state_t*      q;
        residual_state_t* rhs;
        time_state_t      t;
        time_state_t      dt;
        const rhs_calc_t*        rhs_calc;
        const implicit_solve_t*  solver;
        const state_trans_t*     ftrans;
        const inv_state_trans_t* itrans;
        ctrs::array<coeff_t,order+1> diff_coeffs = finite_diff::backward_difference_coeffs_node_based<coeff_t, order>();
        
        detail::cascade_t<var_state_t, order, policy> auxiliary_states;
        
        bdf_t(
            var_state_t& q_in,
            residual_state_t& rhs_in,
            const time_state_t& t0_in,
            const time_state_t& dt_in,
            const rhs_calc_t& rhs_calc_in,
            const implicit_solve_t& solver_in,
            const static_math::int_const_t<order>& order_in = static_math::int_const_t<default_bdf_order>(),
            const state_trans_t&     ftrans = identity_transform,
            const inv_state_trans_t& itrans = identity_transform
        )
        {
            q        = &q_in; //assume filled with initial condition
            rhs      = &rhs_in;
            t        = t0_in;
            dt       = dt_in;
            solver   = &solver_in;
            rhs_calc = &rhs_calc_in;
            solver   = &solver_in;
            ctrs::fill_array(auxiliary_states, (*q));
        }
        
        void advance()
        {
            //the "furthest back" variable state gets set to the residual term            
            auto& grouped_terms = auxiliary_states[0];
            grouped_terms *= diff_coeffs[0];
            for (auto i: range(1,auxiliary_states.size())) {}
            //something
            grouped_terms /= 
            *solver
        }
        
        time_state_t& time(void)     {return  t;}
        var_state_t&  solution(void) {return *q;}
        
        // var_state_t* principal_state;
        // detail::var_state_cascade_t<var_state_t,> aux_states[order];
    };
}