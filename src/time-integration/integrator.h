#pragma once

#include "time-integration/time_integration_utils.h"

namespace spade::time_integration
{
    struct identity_transform_t
    {
        template <typename input_t> void transform_forward(const input_t& input) const {}
        template <typename input_t> void transform_inverse(const input_t& input) const {}
    };
    
    template <
        typename time_value_t,
        typename data_t,
        typename rhs_calc_t,
        typename state_trans_t = identity_transform_t>
    struct integrator_t
    {
        time_axis_t<time_value_t> axis;
        data_t& data;
        const rhs_calc_t& rhs_calc;
        
        integrator_t(
            const time_axis_t<time_value_t>& axis_in,
            data_t& data_in,
            const rhs_calc_t& rhs_calc_in
        )
             : axis{axis_in},
             data{data_in}
             rhs_calc{rhs_calc_in}
        {
            
        }
    };
}