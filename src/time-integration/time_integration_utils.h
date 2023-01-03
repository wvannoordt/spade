#pragma once

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
    
    template <typename>
}