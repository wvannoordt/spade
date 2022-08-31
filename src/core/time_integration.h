#pragma once

#include <concepts>
#include "core/finite_diff.h"
#include "core/iterative_control.h"
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
        time_state_t* t;
        time_state_t  dt;
        var_state_t* q;
        rhs_state_t* k0;
        rhs_state_t* k1;
        rhs_state_t  k1_storage;
        const rhs_calc_t*  rhs_oper;
        const state_trans_t* state_trans;
        const inv_state_trans_t* inv_state_trans;
        rk2(){}
        rk2(
            var_state_t& q_in,
            rhs_state_t& k0_in,
            time_state_t& t_in,
            const time_state_t& dt_in,
            const rhs_calc_t& rhs_oper_in,
            const state_trans_t& trans_in,
            const inv_state_trans_t& inv_trans_in)
        {
            dt = dt_in;
            q  = &q_in;
            t  = &t_in;
            k0 = &k0_in;
            k1_storage = (*k0);
            k1 = &k1_storage;
            rhs_oper = &rhs_oper_in;
            state_trans = &trans_in;
            inv_state_trans = &inv_trans_in;
            print(__FILE__, __LINE__, t);
        }
        
        static std::size_t num_substeps(void) {return 2;}
        
        void advance()
        {
            print(__FILE__, __LINE__, t);
            (*rhs_oper)((*k0), (*q), *t);
            (*k0) *= (0.5*dt);
            (*state_trans)((*q));
            (*q) += (*k0);
            (*inv_state_trans)((*q));
            *t += 0.5*dt;
            (*rhs_oper)((*k1), (*q), *t);
            (*k1) *= (dt);
            (*state_trans)((*q));
            (*q) -= (*k0);
            (*q) += (*k1);
            (*inv_state_trans)((*q));
            *t += 0.5*dt;
        }
        
        time_state_t& time(void)     {return  t;}
        var_state_t&  solution(void) {return *q;}
    };

    namespace detail
    {
        template
        <
            typename coeffs_t,
            typename var_state_t,
            typename time_state_t,
            typename rhs_calc_t,
            typename state_trans_t,
            typename inv_state_trans_t
        >
        struct bdf_inner_rhs_helper
        {
            const coeffs_t* diff_coeffs;
            const var_state_t* grouped_terms;
            time_state_t* t;
            time_state_t dt;
            const rhs_calc_t* rhs_calc;
            const state_trans_t* ftrans;
            const inv_state_trans_t* itrans;
            
            bdf_inner_rhs_helper(
                coeffs_t& diff_coeffs_in,
                var_state_t& grouped_terms_in,
                time_state_t& time_in,
                const time_state_t& dt_in,
                const rhs_calc_t& rhs_calc_in,
                const state_trans_t& ftrans_in,
                const inv_state_trans_t& itrans_in
            )
            {
                diff_coeffs = &diff_coeffs_in;
                grouped_terms = &grouped_terms_in;
                t = &time_in;
                dt = dt_in;
                rhs_calc = &rhs_calc_in;
                ftrans = &ftrans_in;
                itrans = &itrans_in;
            }
            template <typename residual_state_t, typename var_state_in_t>
            void operator()(residual_state_t& rhs, var_state_in_t& q, const time_state_t& pseudotime) const
            {
                //note: using t_n+1 instead of the pseudotime here.
                rhs = (typename coeffs_t::value_type)(0);
                time_state_t cur_time = (*t) + (dt);
                (*rhs_calc)(rhs, q, cur_time);
                rhs *= (dt)/(*diff_coeffs)[(*diff_coeffs).size()-1];
                rhs += (*grouped_terms);
                (*ftrans)(q);
                rhs += q;
                (*itrans)(q);
            }
        };
    }

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
        time_state_t*     t;
        time_state_t      dt;
        const rhs_calc_t*        rhs_calc;
        implicit_solve_t*  solver;
        const state_trans_t*     ftrans;
        const inv_state_trans_t* itrans;
        ctrs::array<coeff_t,order+1> diff_coeffs = finite_diff::backward_difference_coeffs_node_based<coeff_t, order>();
        
        detail::cascade_t<var_state_t, order, policy> auxiliary_states;
        
        bdf_t(){}
        bdf_t(
            var_state_t& q_in,
            residual_state_t& rhs_in,
            time_state_t& t0_in,
            const time_state_t& dt_in,
            const rhs_calc_t& rhs_calc_in,
            implicit_solve_t& solver_in,
            const static_math::int_const_t<order>& order_in = static_math::int_const_t<default_bdf_order>(),
            const state_trans_t&     ftrans_in = identity_transform,
            const inv_state_trans_t& itrans_in = identity_transform
        )
        {
            q        = &q_in; //assume filled with initial condition
            rhs      = &rhs_in;
            t        = &t0_in;
            dt       = dt_in;
            solver   = &solver_in;
            rhs_calc = &rhs_calc_in;
            solver   = &solver_in;
            ftrans   = &ftrans_in;
            itrans   = &itrans_in;
            ctrs::fill_array(auxiliary_states, (*q));
            for (auto i: range(0,auxiliary_states.size())) (*ftrans)(auxiliary_states[i]);
        }
        
        var_state_t& get_grouped_term()
        {
            return auxiliary_states[0];
        }
        
        //We need to compute the first few states with a reduced-order BDF scheme!
        void initialize_states()
        {
            //TODO
        }
        
        void advance()
        {
            //the "furthest back" variable state gets set to the residual term
            auto& grouped_terms = this->get_grouped_term();
            grouped_terms *= diff_coeffs[0]/diff_coeffs[1];
            for (auto i: range(1,auxiliary_states.size()))
            {
                grouped_terms+=auxiliary_states[i];
                grouped_terms*=(diff_coeffs[i]/diff_coeffs[i+1]);
            }
            detail::bdf_inner_rhs_helper inner_rhs(
                diff_coeffs,
                grouped_terms,
                *t,
                dt,
                *rhs_calc,
                *ftrans,
                *itrans
            );
            (*solver)(*rhs, *q, inner_rhs);
            (*ftrans)(*q);
            auxiliary_states.push(*q);
            (*itrans)(*q);
            *t += dt;
        }
        
        time_state_t& time(void)     {return  *t;}
        var_state_t&  solution(void) {return *q;}
        
        // var_state_t* principal_state;
        // detail::var_state_cascade_t<var_state_t,> aux_states[order];
    };
    
    namespace detail
    {
        template
        <
            typename time_integrator_t,
            typename convergence_t,
            typename var_state_t,
            typename rhs_state_t,
            typename rhs_calc_t
        > struct dual_time_inner_update_helper
        {
            time_integrator_t* inner_scheme;
            const convergence_t* convergence_crit;
            dual_time_inner_update_helper(
                time_integrator_t& inner_scheme_in,
                const convergence_t& convergence_crit_in
            )
            {
                inner_scheme = &inner_scheme_in;
                convergence_crit = &convergence_crit_in;
            }
            
            // steady solver
            void operator()(auto& rhs_in, auto& q_in, const auto& rhs_calc)
            {
                using err_t =  typename convergence_t::err_t;
                err_t eps = (1000.0)*convergence_crit->error_tol;
                int num_its =  0;
                // while((eps>convergence_crit->error_tol) && (num_its++ < convergence_crit->max_its))
                while((eps>convergence_crit->error_tol))
                // while(num_its++ < convergence_crit->max_its)
                {
                    inner_scheme->advance();
                    eps = convergence_crit->calc_residual(rhs_in);
                }
            }
        };
    }

    template
    <
        typename var_state_t,
        typename residual_state_t,
        typename time_state_t,
        typename rhs_calc_t,
        typename convergence_t,
        const int bdf_order=default_bdf_order,
        typename state_trans_t     = identity_transform_t,
        typename inv_state_trans_t = identity_transform_t
    > struct dual_time_t
    {
        
        var_state_t* q;
        residual_state_t* rhs;
        time_state_t* t;
        time_state_t  dt;
        const rhs_calc_t* rhs_calc;
        const convergence_t* convergence_crit;
        const state_trans_t*     ftrans;
        const inv_state_trans_t* itrans;
        
        typedef typename detail::get_fundamental_type_if_avail<time_state_t, var_state_t>::type coeff_t;
        typedef detail::bdf_inner_rhs_helper
        <
            ctrs::array<coeff_t,bdf_order+1>,
            var_state_t,
            time_state_t,
            rhs_calc_t,
            state_trans_t,
            inv_state_trans_t
        > inner_rhs_t;
        
        typedef rk2
        <
            var_state_t,
            residual_state_t,
            time_state_t,
            inner_rhs_t,
            state_trans_t,
            inv_state_trans_t
        > inner_scheme_t;
        
        typedef detail::dual_time_inner_update_helper
        <
            inner_scheme_t,
            convergence_t,
            var_state_t,
            residual_state_t,
            rhs_calc_t
        > outer_helper_t;
        
        typedef bdf_t
        <
            var_state_t,
            residual_state_t,
            time_state_t,
            rhs_calc_t,
            outer_helper_t,
            bdf_order,
            state_trans_t,
            inv_state_trans_t
        > outer_scheme_t;
        
        inner_scheme_t inner_scheme;
        
        outer_scheme_t outer_scheme;
        
        dual_time_t(
            var_state_t& q_in,
            residual_state_t& rhs_in,
            time_state_t& t_in,
            const time_state_t& dt_in,
            const rhs_calc_t& rhs_calc_in,
            const convergence_t& convergence_crit_in,
            const static_math::int_const_t<bdf_order>& order_in = static_math::int_const_t<default_bdf_order>(),
            const state_trans_t&     ftrans_in = identity_transform,
            const inv_state_trans_t& itrans_in = identity_transform
        )
        {
            q   = &q_in;
            rhs = &rhs_in;
            t   = &t_in;
            dt  = dt_in;
            rhs_calc = &rhs_calc_in;
            convergence_crit = &convergence_crit_in;
            outer_helper_t inner_solver(inner_scheme, *convergence_crit);
            outer_scheme = outer_scheme_t(*q, *rhs, *t, dt, *rhs_calc, inner_solver, order_in, ftrans_in, itrans_in);
            inner_rhs_t inner_rhs(
                outer_scheme.diff_coeffs,
                outer_scheme.get_grouped_term(),
                *t,
                dt,
                rhs_calc_in,
                ftrans_in,
                itrans_in
            );
            inner_scheme = inner_scheme_t(
                *q,
                *rhs,
                *t,
                dt,
                inner_rhs,
                ftrans_in,
                itrans_in
            );
        }
        
        outer_scheme_t& get_outer_scheme() {return outer_scheme;}
        inner_scheme_t& get_inner_scheme() {return inner_scheme;}
        time_state_t& time(void)     {return *t;}
        var_state_t&  solution(void) {return *q;}
        
        void advance()
        {
            outer_scheme.advance();
        }
    };
}