#pragma once

#include "time-integration/time_integration_utils.h"

namespace spade::time_integration
{
    template <
        typename time_value_t,
        typename scheme_t,
        typename data_t,
        typename rhs_calc_t,
        typename state_trans_t>
    struct integrator_t
    {
        time_axis_t<time_value_t> axis;
        const scheme_t& scheme;
        data_t& data;
        integrator_t(
            const time_axis_t<time_value_t>& axis_in,
            const scheme_t& scheme_in,
            data_t& data_in,
            rhs_calc_t
        )
             : axis{axis_in},
             scheme{scheme_in},
             data{data_in}
        {
            
        }
    };
}