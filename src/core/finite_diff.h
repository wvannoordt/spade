#pragma once

#include "core/config.h"
#include "core/static_math.h"
#include "core/ctrs.h"
#include "core/static_for.h"

namespace spade::finite_diff
{
    template <std::floating_point coeff_t, const std::size_t order>
    static constexpr ctrs::array<coeff_t, order+1> backward_difference_coeffs_node_based()
    {
        ctrs::array<coeff_t, order+1> output;
        coeff_t last = coeff_t(0.0);
        algs::static_for<1,order+1>([&](auto i) -> void
        {
            const coeff_t fact = (coeff_t)static_math::factorial<order>::value/((coeff_t)(static_math::factorial<order-i.value>::value)*((coeff_t)static_math::factorial<i.value>::value));
            output[order-i.value] = ((coeff_t)(1.0)/i.value)*static_math::pow<-1,i.value-1>::value*fact;
            last-=output[order-i.value];
        });
        output[order] = last;
        return output;
    }
    
    template <const int idx, const int order, typename float_t>
    struct centered_finite_diff_t
    {
        static_assert(order == 2*(order/2), "even nominal order required for centered FD coefficients");
        static_assert((idx <= order/2) && (-order/2 <= idx), "index for cenetered FD coefficient out of range [-N/2, N/2].");
        constexpr static int     h     = order/2;
        constexpr static float_t d0    = float_t(static_math::factorial<h>::value)/float_t(static_math::factorial<h+idx>::value);
        constexpr static float_t d1    = float_t(static_math::factorial<h>::value)/float_t(static_math::factorial<h-idx>::value);
        constexpr static float_t value = d0*d1*static_math::pow< -1, 1+idx>::value/float_t(idx);
    };
    
    template <const int order, typename float_t>
    struct centered_finite_diff_t<0, order, float_t>
    {
        constexpr static float_t value = 0.0;
    };
}