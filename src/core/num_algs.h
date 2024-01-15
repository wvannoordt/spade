#pragma once

#include <tuple>
#include <concepts>

#include "core/utils.h"

namespace spade::num_algs
{
    template <typename data_t, std::floating_point float_t>
    struct newton_result_t
    {
        data_t  x;
        int     its;
        float_t eps;
    };
    
    template <std::floating_point data_t>
    struct num_deriv_t
    {
        data_t eps;
    };

    template <std::floating_point data_t>    
    _sp_hybrid inline auto num_deriv(const data_t& eps)
    {
        return num_deriv_t<data_t>{eps};
    }
    
    template <
        typename               data_t,
        std::invocable<data_t> fx_t,
        std::invocable<data_t> dfdx_t,
        std::floating_point    float_t,
        std::invocable<data_t> err_t>
    _sp_hybrid inline auto newton(const data_t& x0, const fx_t& fx, const dfdx_t& dfdx, const int max_its, const float_t& tol, const err_t& err)
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
        return newton_result_t<data_t, float_t>{x, its, eps};
    }
    
    template <
        std::floating_point    data_t,
        std::invocable<data_t> fx_t,
        std::invocable<data_t> dfdx_t,
        std::floating_point    float_t>
    _sp_hybrid inline auto newton(const data_t& x0, const fx_t& fx, const dfdx_t& dfdx, const int max_its, const float_t& tol)
    {
        return newton(x0, fx, dfdx, max_its, tol, [](const data_t& x) { return utils::abs(x);});
    }
    
    template <
        std::floating_point    data_t,
        std::invocable<data_t> fx_t,
        std::floating_point    float_t>
    _sp_hybrid inline auto newton(const data_t& x0, const fx_t& fx, const num_deriv_t<data_t>& nderiv, const int max_its, const float_t& tol)
    {
        return newton(x0, fx, [&](const data_t& x){ return (fx(x+nderiv.eps)-fx(x-nderiv.eps))/(2*nderiv.eps); }, max_its, tol, [](const data_t& x) { return utils::abs(x);});
    }
}