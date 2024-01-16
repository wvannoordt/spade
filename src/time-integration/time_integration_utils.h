#pragma once

#include "core/ctrs.h"

#include "time-integration/explicit.h"
#include "time-integration/implicit.h"

namespace spade::time_integration
{
    template <typename time_value_t> struct time_axis_t
    {
        using value_type = time_value_t;
        time_value_t t, dt;
        time_axis_t(const time_value_t& t0_in, const time_value_t& dt_in): t{t0_in}, dt{dt_in}
        {}
        const time_value_t& timestep() const {return dt;}
        time_value_t& timestep() {return dt;}
        const time_value_t& time() const {return t;}
        time_value_t& time() {return t;}
    };
    
    template <typename var_state_t, typename rhs_state_t, typename scheme_t>
    struct integrator_data_t
    {
        ctrs::array<var_state_t, scheme_t::var_size()> solution_data;
        ctrs::array<rhs_state_t, scheme_t::rhs_size()> residual_data;
        const scheme_t scheme_data;
        
        using index_type = int;
        using solution_type = var_state_t;
        using residual_type = rhs_state_t;
        using scheme_type   = scheme_t;
        
        //For memory-optimization / external initialization
        integrator_data_t(const scheme_t& scheme_in) : scheme_data{scheme_in} {}
        
        integrator_data_t(const var_state_t& q, const rhs_state_t& r, const scheme_t& scheme_in) : scheme_data{scheme_in}
        {
            //need to improve this
            solution_data = q;
            residual_data = r;
        }
        
        integrator_data_t(var_state_t&& q, rhs_state_t&& r, const scheme_t& scheme_in) : scheme_data{scheme_in}
        {
            solution_data[0] = std::move(q);
            residual_data[0] = std::move(r);
            for (int i = 1; i < solution_data.size(); ++i) solution_data[i] = solution_data[0];
            for (int i = 1; i < residual_data.size(); ++i) residual_data[i] = residual_data[0];
        }
        
        const rhs_state_t& residual(const index_type i = 0) const {return residual_data[i];}
        rhs_state_t& residual(const index_type i = 0) {return residual_data[i];}
        const var_state_t& solution(const index_type i = 0) const {return solution_data[i];}
        var_state_t& solution(const index_type i = 0) {return solution_data[i];}
        const scheme_t& scheme() const {return scheme_data;}
    };
    
    //Need to modify this to allow for large copyable values
    template <typename... integrator_datas_t>
    auto couple(integrator_datas_t&... datas)
    {
        return std::forward_as_tuple(datas...);
    }
}