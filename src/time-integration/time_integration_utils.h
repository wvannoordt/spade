#pragma once

#include "core/ctrs.h"

#include "time-integration/explicit.h"
#include "time-integration/implicit.h"

namespace spade::time_integration
{
    template <typename time_value_t> struct time_axis_t
    {
        time_value_t t0, dt;
        time_axis_t(const time_value_t& t0_in, const time_value_t& dt_in): t0{t0_in}, dt{dt_in}
        {}
    };
    
    template <typename var_state_t, typename rhs_state_t, typename scheme_t>
    struct integrator_data_t
    {
        ctrs::array<var_state_t, scheme_t::var_size()> solution_data;
        ctrs::array<rhs_state_t, scheme_t::rhs_size()> residual_data;
        const scheme_t& scheme_data;
        
        using index_type = int;
        
        integrator_data_t(const var_state_t& q, const rhs_state_t& r, const scheme_t& scheme_in) : scheme_data{scheme_in}
        {
            solution_data = q;
            residual_data = r;
        }
        
        const rhs_state_t& residual(const index_type i = 0) const {return residual_data[i];}
        rhs_state_t& residual(const index_type i = 0) {return residual_data[i];}
        const var_state_t& solution(const index_type i = 0) const {return solution_data[i];}
        var_state_t& solution(const index_type i = 0) {return solution_data[i];}
        const scheme_t& scheme() const {return scheme_data;}
    };
}