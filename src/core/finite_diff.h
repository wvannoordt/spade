#pragma once

#include "core/config.h"
#include "core/static_math.h"
#include "core/ctrs.h"

namespace spade::finite_diff
{
    template <std::floating_point coeff_t, const std::size_t order>
    static constexpr ctrs::array<coeff_t, order+1> backward_difference_coeffs_node_based()
    {
        ctrs::array<coeff_t, order+1> output;
        coeff_t last = coeff_t(0.0);
        static_for<1,order+1>([&](auto i) -> void
        {
            const coeff_t fact = (coeff_t)static_math::factorial<order>::value/((coeff_t)(static_math::factorial<order-i.value>::value)*((coeff_t)static_math::factorial<i.value>::value));
            output[order-i.value] = ((coeff_t)(1.0)/i.value)*static_math::pow<-1,i.value-1>::value*fact;
            last-=output[order-i.value];
        });
        output[order] = last;
        return output;
    }
}