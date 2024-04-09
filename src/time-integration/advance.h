#pragma once

#include <concepts>
#include <ratio>

#include "time-integration/time_integration_utils.h"
#include "time-integration/explicit.h"
#include "time-integration/implicit.h"

namespace spade::time_integration
{
    namespace detail
    {
        using ratio_integral_type = decltype(std::deci::num);
        
        template <typename val_t> struct nonzero_t
        {
            constexpr static bool value = false;
        };
        
        template <const ratio_integral_type num, const ratio_integral_type denum> struct nonzero_t<std::ratio<num, denum>>
        {
            constexpr static bool value = (num!=0);
        };
        
        template <const ratio_integral_type num> struct nonzero_t<std::ratio<num>>
        {
            constexpr static bool value = (num!=0);
        };
        
        template <typename numeric_t, typename type_t> struct coeff_value_t
        {
            //invalid
        };
        
        template <typename numeric_t, const ratio_integral_type num, const ratio_integral_type denum> struct coeff_value_t<numeric_t, std::ratio<num, denum>>
        {
            constexpr static numeric_t value() {return numeric_t(num)/numeric_t(denum);}
        };
        
        template <typename numeric_t, const ratio_integral_type num> struct coeff_value_t<numeric_t, std::ratio<num>>
        {
            constexpr static numeric_t value() {return numeric_t(num);}
        };
        
        template <typename idx_t> const int idx_val = utils::remove_all<idx_t>::type::value;
        
        template <typename a_t, typename b_t> struct ratio_diff_t;
        
        template <ratio_integral_type n0, ratio_integral_type d0, ratio_integral_type n1, ratio_integral_type d1>
        struct ratio_diff_t<std::ratio<n0, d0>, std::ratio<n1, d1>>
        {
            using type = std::ratio<d1*n0 - d0*n1, d0*d1>;
        };
        
        template <typename prev_row_t, typename curr_row_t, typename dt_t, typename sol_t, typename gas_t, typename state_t, const std::size_t nstage, typename rhs_t>
        static void transform_advance_to(
            const prev_row_t&,
            const curr_row_t&,
            const dt_t& dt,
            sol_t& q,
            const fluid_state::state_transform_t<gas_t, state_t>& trans,
            const ctrs::array<rhs_t, nstage>& resids)
        {
            using numeric_type = dt_t;
            using rhs_img_type = typename rhs_t::const_image_type;
            ctrs::array<rhs_img_type, nstage> rhs_images;
            for (int j = 0; j < nstage; ++j) rhs_images[j] = resids[j].image();
            auto q_img = q.image();
            
            using alias_type = typename sol_t::alias_type;
            using resid_type = typename rhs_t::alias_type;
            const gas_t gas_model = trans.gas;
            
            algs::transform_inplace(q, [=] _sp_hybrid (const grid::cell_idx_t& ii, const alias_type& q_orig)
            {
                state_t cons;
                fluid_state::convert_state(q_orig, cons, gas_model);
                
                constexpr int nupdate = curr_row_t::length();
                algs::static_for<0, nupdate>([&](const auto& idx)
                {
                    constexpr int i = idx.value;
                    const auto rhs_img = rhs_images[i];
                    using prev_coeff_t = typename prev_row_t::elem_t<i>;
                    using next_coeff_t = typename curr_row_t::elem_t<i>;
                    using diff_coeff_t = typename ratio_diff_t<next_coeff_t, prev_coeff_t>::type;
                    constexpr bool nonzero_coeff = nonzero_t<diff_coeff_t>::value;
                    const numeric_type dt_tmp = dt;
                    if constexpr (nonzero_coeff)
                    {
                        numeric_type time_coeff = detail::coeff_value_t<numeric_type, diff_coeff_t>::value()*dt_tmp;
                        const resid_type rhs_elem = rhs_img.get_elem(ii);
                        cons += time_coeff*rhs_elem;
                    }
                });
                
                alias_type new_q;
                fluid_state::convert_state(cons, new_q, gas_model);
                return new_q;
            });
        }
    }
    
    
    //explicit RK scheme (TODO: implement concept for this)
    template <typename axis_t, typename data_t, typename scheme_t, typename rhs_t, typename boundary_t, typename trans_t>
    requires (data_t::scheme_type::is_rk_specialization)
    void integrate_advance(axis_t& axis, data_t& data, const scheme_t& scheme, const rhs_t& rhs, const boundary_t& boundary, const trans_t& trans)
    {
        //Note that we assume that the boundary has been set before we begin integration.
        //Note also that the boundary pertains to the untransformed state
        
        constexpr bool has_second_variable_copy = data.solution_data.size() > 1;

        //Number of RHS evaluations
        const int num_stages = scheme_t::table_type::rows();
        
        //Get the basic numeric type from the solution
        using numeric_type = typename axis_t::value_type;
        
        //Timestep
        const auto& dt = axis.timestep();
        //Computing the residuals
        algs::static_for<0, num_stages>([&](const auto& ii) -> void
        {
            const int i = ii.value;
            //We won't update this copy of the solution during the RHS evaluations
            // const auto& sol_base = data.solution(0);
            
            //The solution used when we evaluate the RHS
            constexpr int sol_idx = has_second_variable_copy ? 1 : 0;
            auto& sol = data.solution(sol_idx);
            if constexpr (has_second_variable_copy)
            {
                sol = data.solution(0);
            }
            // sol = sol_base;//expensive!
            
            //solution augmentation loop
            //Begin by applying the forward transform
            if constexpr (i != 0) trans.transform_forward(sol);
            algs::static_for<0, i>([&](const auto& jj) -> void
            {
                const int j = jj.value;
                using table_t = scheme_t::table_type;
                //Note that this is a temporary workaround owing to a bug in GCC 10.2
                using coeff_value_t_1 = typename table_t::elem_t<detail::idx_val<decltype(ii)>>;
                using coeff_value_t = typename coeff_value_t_1::elem_t<detail::idx_val<decltype(jj)>>;
                
                //Using GCC 11+, below is valid
                //using coeff_value_t = typename table_t::elem_t<i>::elem_t<j>;
                
                if constexpr (detail::nonzero_t<coeff_value_t>::value)
                {
                    constexpr numeric_type coeff = detail::coeff_value_t<numeric_type, coeff_value_t>::value(); //get the coefficient
                    auto& resid_j = data.residual(j); //"k_j"
                    resid_j *= (dt*coeff); //This requires that "dt" is cheap to multiply, otherwise there is a more efficient way of doing this!
                    sol += resid_j;
                    resid_j *= numeric_type(1.0)/(dt*coeff);
                }
            });
            if constexpr (i != 0) trans.transform_inverse(sol);
            //By now, the solution has been augmented to accommodate the
            //evaluation of the ith residual
            auto& cur_resid = data.residual(i);
            using dt_coeff_value_t = typename scheme_t::dt_type::elem_t<i>;
            constexpr numeric_type time_coeff = detail::coeff_value_t<numeric_type, dt_coeff_value_t>::value();
            
            //Evaluate the residual at t + c*dt
            axis.time() += time_coeff*dt;
            if constexpr (i > 0) boundary(sol, axis.time());
            rhs(cur_resid, sol, axis.time());
            axis.time() -= time_coeff*dt;
            
            //solution augmentation loop
            //Begin by applying the forward transform
            if constexpr (!has_second_variable_copy)
            {
                trans.transform_forward(sol);
                algs::static_for<0, i>([&](const auto& jj) -> void
                {
                    const int j = jj.value;
                    using table_t = scheme_t::table_type;
                    //Note that this is a temporary workaround owing to a bug in GCC 10.2
                    using coeff_value_t_1 = typename table_t::elem_t<detail::idx_val<decltype(ii)>>;
                    using coeff_value_t = typename coeff_value_t_1::elem_t<detail::idx_val<decltype(jj)>>;
                    
                    //Using GCC 11+, below is valid
                    //using coeff_value_t = typename table_t::elem_t<i>::elem_t<j>;
                    
                    if constexpr (detail::nonzero_t<coeff_value_t>::value)
                    {
                        constexpr numeric_type coeff = detail::coeff_value_t<numeric_type, coeff_value_t>::value(); //get the coefficient
                        auto& resid_j = data.residual(j); //"k_j"
                        resid_j *= (dt*coeff); //This requires that "dt" is cheap to multiply, otherwise there is a more efficient way of doing this!
                        sol -= resid_j;
                        resid_j *= numeric_type(1.0)/(dt*coeff);
                    }
                });
                trans.transform_inverse(sol);
            }
        });
        auto& new_solution = data.solution(0); //solution is updated in place
        
        //Residual accumulation loop
        //Update the transformed solution
        trans.transform_forward(new_solution);
        algs::static_for<0, num_stages>([&](const auto& ii) -> void
        {
            const int i = ii.value;
            using accum_coeff_value_t = typename scheme_t::accum_type::elem_t<i>;
            
            //We skip over any zero coefficients
            if constexpr (detail::nonzero_t<accum_coeff_value_t>::value)
            {
                constexpr numeric_type coeff = detail::coeff_value_t<numeric_type, accum_coeff_value_t>::value();
                auto& resid = data.residual(i);
                //Multiply the residual by the accumulation coefficient
                resid *= (coeff*dt);
                
                new_solution += resid;
                //Divide to avoid modifying the residual unduly (may be unnecessary!)
                resid *= numeric_type(1.0)/(coeff*dt);
            }
        });
        trans.transform_inverse(new_solution);
        axis.time() += dt;
        boundary(new_solution, axis.time());
    }
    
    
    //Experimental version
    template <typename axis_t, typename data_t, typename scheme_t, typename rhs_t, typename boundary_t, typename state_t, typename gas_t>
    requires (data_t::scheme_type::is_rk_specialization)
    void integrate_advance(axis_t& axis, data_t& data, const scheme_t& scheme, const rhs_t& rhs, const boundary_t& boundary, const fluid_state::state_transform_t<gas_t, state_t>& trans)
    {
        //Note that we assume that the boundary has been set before we begin integration.
        auto& q = data.solution(0);
        
        //Number of RHS evaluations
        const int num_stages = scheme_t::table_type::rows();
        
        //Get the basic numeric type from the solution
        using numeric_type = typename axis_t::value_type;
        
        //Timestep
        const auto& dt = axis.timestep();
        
        auto& rhs0 = data.residual(0);
        using dt_coeff_value_t = typename scheme_t::dt_type::elem_t<0>;
        constexpr numeric_type time_coeff = detail::coeff_value_t<numeric_type, dt_coeff_value_t>::value();
        axis.time() += time_coeff*dt;
        rhs(rhs0, q, axis.time());
        axis.time() -= time_coeff*dt;
        
        algs::static_for<1, num_stages>([&](const auto& i_substep)
        {
            constexpr int i = i_substep.value;
            using prev_row_t = typename scheme_t::table_type::elem_t<i - 1>;
            using curr_row_t = typename scheme_t::table_type::elem_t<i>;
            detail::transform_advance_to(prev_row_t(), curr_row_t(), dt, q, trans, data.residual_data);
            
            using dt_coeff_value_t = typename scheme_t::dt_type::elem_t<i>;
            constexpr numeric_type time_coeff = detail::coeff_value_t<numeric_type, dt_coeff_value_t>::value();
            
            axis.time() += time_coeff*dt;
            boundary(q, axis.time());
            rhs(data.residual(i), q, axis.time());
            axis.time() -= time_coeff*dt;
        });
        
        using last_row_t = typename scheme_t::table_type::elem_t<num_stages - 1>;
        using resd_row_t = typename scheme_t::accum_type;
        detail::transform_advance_to(last_row_t(), resd_row_t(), dt, q, trans, data.residual_data);
        
        //After it all, we set boundary conditions
        axis.time() += dt;
        boundary(q, axis.time());
    }
    
    
    template <typename axis_t, typename data_t, typename scheme_t, typename rhs_t, typename boundary_t, typename trans_t>
    requires (data_t::scheme_type::is_crank_nichol_specialization)
    void integrate_advance(axis_t& axis, data_t& data, const scheme_t& scheme, const rhs_t& rhs, const boundary_t& boundary, const trans_t& trans)
    {
        //Note that we assume that the boundary has been set before we begin integration.
        //Note also that the boundary pertains to the untransformed state

        //Number of RHS evaluations
        auto& q   = data.solution(0);
        auto& del = data.residual(0);
        auto& res = data.residual(1);

        using coeff_t          = typename data_t::scheme_type::coeff_type;
        const coeff_t beta     = (coeff_t(1.0) - data.scheme_data.coeff);
        const coeff_t alpha    = data.scheme_data.coeff;
        
        const auto& dt = axis.timestep();
        const coeff_t gamma = coeff_t(0.5)*dt;
        
        rhs(del, q, axis.time());
        del *= gamma;
        trans.transform_forward(q);
        del += q;
        trans.transform_inverse(q);
        
        bool converged = false;
        int num_its = 0;
        axis.time() += dt;
        
        coeff_t err = coeff_t(1e10);
        
        while (!converged)
        {                                 //  q         del                               res
            rhs(res, q, axis.time());     //  qi        w0 + gam*R0                       Ri
            res *= gamma;                 //                                              gam*Ri
            del += res;                   //            w0 + gam*R0 + gam*Ri
            trans.transform_forward(q);   //  wi
            del -= q;                     //            w0 - wi + gam*R0 + gam*Ri
            del *= beta;                  //            bet*(w0 - wi + gam*R0 + gam*Ri)
            res -= q;                     //                                              gam*Ri - wi
            q   += del;                   //  w{i+1}
            err  = scheme.err_func(del);  // (compute |w{i+1} - wi|)
            trans.transform_inverse(q);   //  q{i+1}
            del *= coeff_t(1.0)/beta;     //            w0 - wi + gam*R0 + gam*Ri
            del -= res;                   //            w0 + gam*R0
            
            boundary(q, axis.time());
            
            num_its++;
            converged = ((err < scheme.crit.tolerance()) || (num_its > scheme.crit.max_its()));
        }
    }
}