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
        template <typename T> concept has_fundamental_type = requires(T t) {T::fundamental_type;} && std::floating_point<typename T::fundamental_type>;
        template <typename query_t> struct get_numeric_type{};
        template <std::floating_point query_t>  struct get_numeric_type<query_t> { using type = query_t; };
        template <has_fundamental_type query_t> struct get_numeric_type<query_t> { using type = typename query_t::fundamental_type; };
        
        template <typename val_t> struct nonzero_t
        {
            constexpr static bool value = false;
        };
        
        template <const ratio_integral_type num, const ratio_integral_type denum> struct nonzero_t<std::ratio<num, denum>>
        {
            constexpr static bool value = (num==0);
        };
        
        template <const ratio_integral_type num> struct nonzero_t<std::ratio<num>>
        {
            constexpr static bool value = (num==0);
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
    }
    
    //explicit RK scheme (TODO: implement concept for this)
    template <typename axis_t, typename data_t, typename scheme_t, typename rhs_t, typename trans_t>
    void integrate_advance(axis_t& axis, data_t& data, const scheme_t& scheme, const rhs_t& rhs, const trans_t& trans)
    {
        const int num_stages = scheme_t::table_type::rows();
        using numeric_type = detail::get_numeric_type<typename data_t::solution_type>::type;
        const auto& dt = axis.timestep();
        algs::static_for<0, num_stages>([&](const auto& ii) -> void
        {
            const int i = ii.value;
            const auto& sol_base = data.solution(0); //The first solution is the solution at the previous timestep
            
            auto& sol = data.solution(1); //all RK schemes only need a single solution array.
            sol = sol_base; //expensive!
            
            //solution augmentation loop
            algs::static_for<0, i>([&](const auto& jj) -> void
            {
                const int j = jj.value;
                using table_t = scheme_t::table_type;
                using coeff_value_t = typename table_t::elem_t<i>::elem_t<j>;
                if constexpr (detail::nonzero_t<coeff_value_t>::value)
                {
                    constexpr numeric_type coeff = detail::coeff_value_t<numeric_type, coeff_value_t>::value(); //get the coefficient
                    auto& resid_j = data.residual(j); //"k_j"
                    print(coeff);
                    resid_j *= (dt*coeff); //This requires that "dt" is cheap to multiply, otherwise there is a more efficient way of doing this!
                    sol += resid_j;
                    resid_j /= (dt*coeff);
                }
            });
            
            //by now, the solution has been augmented to accommodate the
            //evaluation of the ith residual
            auto& cur_resid = data.residual(i);
            using dt_coeff_value_t = typename scheme_t::dt_type::elem_t<i>;
            constexpr numeric_type time_coeff = detail::coeff_value_t<numeric_type, dt_coeff_value_t>::value();
            print(time_coeff);
            axis.time() += time_coeff*dt;
            rhs(cur_resid, sol, axis.time());
            axis.time() -= time_coeff*dt;
        });
        print(data.residual(0), data.residual(1), data.residual(2), data.residual(3));
        
        //residual accumulation
        auto& new_solution = data.solution(0);
        algs::static_for<0, num_stages>([&](const auto& ii) -> void
        {
            const int i = ii.value;
            const numeric_type coeff = 0.0; //get the accumulation coefficient
            using accum_coeff_value_t = typename scheme_t::accum_type::elem_t<i>;
            if constexpr (detail::nonzero_t<accum_coeff_value_t>::value)
            {
                constexpr numeric_type coeff = detail::coeff_value_t<numeric_type, accum_coeff_value_t>::value();
                print(coeff);
                auto& resid = data.residual(i);
                resid *= (coeff*dt);
                new_solution += resid;
                resid /= (coeff*dt);
            }
        });
        
        axis.time() += dt;
    }
}