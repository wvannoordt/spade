#pragma once

#include <tuple>
#include <concepts>

#include "core/utils.h"

namespace spade::num_algs
{
    template <
        typename               data_t,
        std::invocable<data_t> fx_t,
        std::invocable<data_t> dfdx_t,
        std::floating_point    float_t,
        std::invocable<data_t> err_t>
    auto newton(const data_t& x0, const fx_t& fx, const dfdx_t& dfdx, const int max_its, const float_t& tol, const err_t& err)
    {
        data_t  x   = x0;
        auto    rhs = fx(x);
        float_t eps = err(rhs);
        int     its = 0;
        while((eps > tol) && (its++ < max_its))
        {
            x  -= rhs/dfdx(x);
            rhs = fx(x);
            eps = err(rhs);
        }
        return std::make_tuple(x, its, eps);
    }
    
    template <
        std::floating_point    data_t,
        std::invocable<data_t> fx_t,
        std::invocable<data_t> dfdx_t,
        std::floating_point    float_t>
    auto newton(const data_t& x0, const fx_t& fx, const dfdx_t& dfdx, const int max_its, const float_t& tol)
    {
        return newton(x0, fx, dfdx, max_its, tol, [](const data_t& x) { return utils::abs(x);});
    }
}