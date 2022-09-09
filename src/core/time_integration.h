#pragma once

#include <concepts>
#include "core/finite_diff.h"
#include "core/iterative_control.h"
#include "core/composite_transform.h"
#include "core/ctrs.h"

namespace spade::time_integration
{
    enum var_cascade_policy
    {
        var_cascade_assign,
        var_cascade_inplace
    };
    
    static constexpr var_cascade_policy default_policy = var_cascade_inplace;
    
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
        template <typename input_t> void transform_forward(const input_t& input) const {}
        template <typename input_t> void transform_inverse(const input_t& input) const {}
    };
    
    static identity_transform_t identity_transform;
    
    template <
        typename var_state_t,
        typename rhs_state_t,
        typename time_state_t,
        typename rhs_calc_t,
        typename state_trans_t
        >
    struct rk2
    {
        time_state_t* t;
        time_state_t  dt;
        var_state_t* q;
        rhs_state_t* k0;
        rhs_state_t* k1;
        rhs_state_t  k1_storage;
        rhs_calc_t*  rhs_oper;
        const state_trans_t* state_trans;
        rk2(){}
        rk2(
            var_state_t& q_in,
            rhs_state_t& k0_in,
            time_state_t& t_in,
            const time_state_t& dt_in,
            rhs_calc_t& rhs_oper_in,
            const state_trans_t& trans_in)
        {
            dt = dt_in;
            q  = &q_in;
            t  = &t_in;
            k0 = &k0_in;
            k1_storage = (*k0);
            k1 = &k1_storage;
            rhs_oper = &rhs_oper_in;
            state_trans = &trans_in;
        }
        
        static std::size_t num_substeps(void) {return 2;}
        
        void advance()
        {
            (*rhs_oper)((*k0), (*q), *t);
            (*k0) *= (0.5*dt);
            state_trans->transform_forward(*q);
            (*q) += (*k0);
            state_trans->transform_inverse(*q);
            *t += 0.5*dt;
            (*rhs_oper)((*k1), (*q), *t);
            (*k1) *= (dt);
            state_trans->transform_forward(*q);
            (*q) -= (*k0);
            (*q) += (*k1);
            state_trans->transform_inverse(*q);
            *t += 0.5*dt;
            
            //This is done so that the residual is physically scaled for external use.
            (*k1) /= (dt);
        }
        
        time_state_t& time(void)     {return *t;}
        var_state_t&  solution(void) {return *q;}
        rhs_state_t&  residual(void) {return *k1;}
    };

    namespace detail
    {
        template
        <
            typename coeffs_t,
            typename var_state_t,
            typename residual_state_t,
            typename time_state_t,
            typename rhs_calc_t,
            typename state_trans_t
        >
        struct bdf_inner_rhs_helper
        {
            const coeffs_t* diff_coeffs;
            const var_state_t* grouped_terms;
            time_state_t* t;
            time_state_t dt;
            rhs_calc_t* rhs_calc;
            const state_trans_t* trans;
            bdf_inner_rhs_helper(){}
            bdf_inner_rhs_helper(
                coeffs_t& diff_coeffs_in,
                var_state_t& grouped_terms_in,
                residual_state_t& rhs,
                time_state_t& time_in,
                const time_state_t& dt_in,
                rhs_calc_t& rhs_calc_in,
                const state_trans_t& trans_in
            )
            {
                diff_coeffs = &diff_coeffs_in;
                set_grouped_term(grouped_terms_in);
                t = &time_in;
                dt = dt_in;
                rhs_calc = &rhs_calc_in;
                trans = &trans_in;
            }
            
            void set_grouped_term(const var_state_t& g) {grouped_terms = &g;}
            
            void operator()(residual_state_t& rhs, var_state_t& q, const time_state_t& pseudotime)
            {
                // note: using t_n+1 instead of the pseudotime here.
                rhs = (typename coeffs_t::value_type)(0);
                time_state_t cur_time = (*t) + (dt);
                (*rhs_calc)(rhs, q, cur_time);
                rhs *= (-dt)/(*diff_coeffs)[(*diff_coeffs).size()-1];
                rhs -= (*grouped_terms);
                trans->transform_forward(q);
                rhs -= q;
                trans->transform_inverse(q);
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
        const var_cascade_policy policy=default_policy
    > struct bdf_t
    {
        static_assert(order <= 6, "BDF non-zero-stable for order larger than 6!");
        typedef typename detail::get_fundamental_type_if_avail<time_state_t, var_state_t>::type coeff_t;
        
        var_state_t*      q;
        residual_state_t* rhs;
        time_state_t*     t;
        time_state_t      dt;
        rhs_calc_t*  rhs_calc;
        implicit_solve_t*  solver;
        const state_trans_t* trans;
        ctrs::array<coeff_t,order+1> diff_coeffs = finite_diff::backward_difference_coeffs_node_based<coeff_t, order>();
        
        detail::cascade_t<var_state_t, order, policy> auxiliary_states;
        
        bdf_t(){}
        bdf_t(
            var_state_t& q_in,
            residual_state_t& rhs_in,
            time_state_t& t0_in,
            const time_state_t& dt_in,
            rhs_calc_t& rhs_calc_in,
            implicit_solve_t& solver_in,
            const static_math::int_const_t<order>& order_in = static_math::int_const_t<default_bdf_order>(),
            const state_trans_t&       trans_in = identity_transform
        )
        {
            q        = &q_in; //assume filled with initial condition
            rhs      = &rhs_in;
            t        = &t0_in;
            dt       = dt_in;
            solver   = &solver_in;
            rhs_calc = &rhs_calc_in;
            solver   = &solver_in;
            trans   = &trans_in;
            ctrs::fill_array(auxiliary_states, (*q));
            for (auto i: range(0,auxiliary_states.size())) trans->transform_forward(auxiliary_states[i]);
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
                *rhs,
                *t,
                dt,
                *rhs_calc,
                *trans
            );
            (*solver)(*rhs, *q, inner_rhs);
            trans->transform_forward(*q);
            auxiliary_states.push(*q);
            trans->transform_inverse(*q);
            *t += dt;
        }
        
        time_state_t& time(void)     {return *t;}
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
            convergence_t* convergence_crit;
            dual_time_inner_update_helper(){}
            dual_time_inner_update_helper(
                time_integrator_t& inner_scheme_in,
                convergence_t& convergence_crit_in
            )
            {
                inner_scheme = &inner_scheme_in;
                convergence_crit = &convergence_crit_in;
            }
            
            // steady solver
            void operator()(rhs_state_t& rhs_in, var_state_t& q_in, rhs_calc_t& rhs_calc)
            {
                using err_t =  typename convergence_t::err_t;
                err_t eps = (1000.0)*convergence_crit->error_tol;
                int num_its =  0;
                while((eps>convergence_crit->error_tol) && (num_its++ < convergence_crit->max_its))                
                {
                    inner_scheme->advance();
                    eps = convergence_crit->calc_residual(inner_scheme->residual());
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
        typename state_trans_t = identity_transform_t,
        typename conditioner_t = identity_transform_t
    > struct dual_time_t
    {
        var_state_t* q;
        residual_state_t* rhs;
        residual_state_t  rhs_inner_val;
        time_state_t* t;
        time_state_t  tau;
        time_state_t  dt;
        time_state_t  dt_inner;
        rhs_calc_t* rhs_calc;
        convergence_t* convergence_crit;
        const state_trans_t* trans;
        const conditioner_t* conditioner;
        
        using composite_transform_type = algs::composite_transform_t<state_trans_t, conditioner_t>;
        composite_transform_type composite_transform;
        
        typedef typename detail::get_fundamental_type_if_avail<time_state_t, var_state_t>::type coeff_t;
        typedef detail::bdf_inner_rhs_helper
        <
            ctrs::array<coeff_t,bdf_order+1>,
            var_state_t,
            residual_state_t,
            time_state_t,
            rhs_calc_t,
            state_trans_t//composite_transform_type
        > inner_rhs_t;
        
        typedef rk2
        <
            var_state_t,
            residual_state_t,
            time_state_t,
            inner_rhs_t,
            composite_transform_type
        > inner_scheme_t;
        
        typedef detail::dual_time_inner_update_helper
        <
            inner_scheme_t,
            convergence_t,
            var_state_t,
            residual_state_t,
            inner_rhs_t
        > outer_helper_t;
        
        typedef bdf_t
        <
            var_state_t,
            residual_state_t,
            time_state_t,
            rhs_calc_t,
            outer_helper_t,
            bdf_order,
            state_trans_t
        > outer_scheme_t;
        
        inner_scheme_t inner_scheme;
        
        outer_scheme_t outer_scheme;
        
        outer_helper_t inner_solver;
        
        inner_rhs_t inner_rhs;
        
        dual_time_t(
            var_state_t& q_in,
            residual_state_t& rhs_in,
            time_state_t& t_in,
            const time_state_t& dt_in,
            const time_state_t& dt_inner_in,
            rhs_calc_t& rhs_calc_in,
            convergence_t& convergence_crit_in,
            const static_math::int_const_t<bdf_order>& order_in = static_math::int_const_t<default_bdf_order>(),
            const state_trans_t& trans_in = identity_transform,
            const conditioner_t& conditioner_in = identity_transform
        )
        {
            q   = &q_in;
            rhs = &rhs_in;
            rhs_inner_val = rhs_in;
            t   = &t_in;
            dt  = dt_in;
            tau = time_state_t(0.0);
            dt_inner  = dt_inner_in;
            rhs_calc = &rhs_calc_in;
            convergence_crit = &convergence_crit_in;
            inner_solver = outer_helper_t(inner_scheme, *convergence_crit);
            outer_scheme = outer_scheme_t(*q, *rhs, *t, dt, *rhs_calc, inner_solver, order_in, trans_in);
            composite_transform.list.set_items(trans_in, conditioner_in);
            inner_rhs = inner_rhs_t(
                outer_scheme.diff_coeffs,
                outer_scheme.get_grouped_term(),
                rhs_inner_val,
                *t,
                dt,
                rhs_calc_in,
                trans_in
            );
            inner_scheme = inner_scheme_t(
                *q,
                rhs_inner_val,
                tau,
                dt_inner,
                inner_rhs,
                composite_transform
            );
            
            //get rid of this crap
            inner_scheme.k1 = &(inner_scheme.k1_storage);
        }
        
        outer_scheme_t& get_outer_scheme() {return outer_scheme;}
        inner_scheme_t& get_inner_scheme() {return inner_scheme;}
        time_state_t& time(void)     {return *t;}
        var_state_t&  solution(void) {return *q;}
        
        void advance()
        {
            inner_rhs.set_grouped_term(outer_scheme.get_grouped_term());
            outer_scheme.advance();
        }
    };
}