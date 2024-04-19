#pragma once

#include "time-integration/time_integration_utils.h"
#include "time-integration/advance.h"
#include "time-integration/advance_limiter.h"

namespace spade::time_integration
{
    struct identity_transform_t
    {
        template <typename input_t> void transform_forward(const input_t& input) const {}
        template <typename input_t> void transform_inverse(const input_t& input) const {}
    };
    
    static identity_transform_t identity_transform;

    static struct no_boundary_t { void operator () (auto& q, const auto& t) const {}} no_boundary;

	static struct no_limiter_t { void operator () (auto& dsol, const auto& sol) const {}} no_limiter;
    
    template <
        typename time_value_t,
        typename scheme_t,
        typename data_t,
        typename rhs_calc_t,
        typename boundary_cond_t = no_boundary_t,
        typename state_trans_t = identity_transform_t,
	    typename advance_limiter_t = no_limiter_t>
    struct integrator_t
    {
        time_axis_t<time_value_t>   axis;
        data_t&                     data;
        scheme_t                    scheme;
        const rhs_calc_t&           rhs_calc;
        const boundary_cond_t&      boundary_cond;
        const state_trans_t&        trans;
        const advance_limiter_t&    advance_lim;
		
        integrator_t(
            const time_axis_t<time_value_t>& axis_in,
            const scheme_t& scheme_in,
            data_t& data_in,
            const rhs_calc_t& rhs_calc_in,
            const boundary_cond_t& boundary_cond_in = no_boundary,
            const state_trans_t& trans_in = identity_transform,
			const advance_limiter_t& advance_lim_in = no_limiter
        )
             : axis        {axis_in},
             data          {data_in},
             scheme        {scheme_in},
             rhs_calc      {rhs_calc_in},
             boundary_cond {boundary_cond_in},
     		 trans         {trans_in},
 		     advance_lim   {advance_lim_in}  
        {
            
        }
        
        void advance()
        {
            integrate_advance(axis, data, scheme, rhs_calc, boundary_cond, trans, advance_lim);
        }
        
        data_t::solution_type& solution() {return data.solution(0);}
        const time_value_t& time() {return axis.time();}
    };
}
