#pragma once

namespace cvdf::algs
{
    template <class ltype, class rtype> auto sum(const ltype& l, const rtype& r)
    {
        return l + r;
    }
    
    template <class invokable_t, class param_t>
    static auto reduce_over_params(const invokable_t& func, const param_t& param)
    {
        return func(param);
    }
    
    template <class invokable_t, class reduction_t, class param_t, class... params_t>
    static auto reduce_over_params(const invokable_t& func, const reduction_t& reduce_op, const param_t& param, params_t... params)
    {
        return reduce_op(func(param),reduce_over_params(func, params...));
    }
    template <class invokable_t, class param_t, class... params_t>
    static auto reduce_over_params(const invokable_t& func, const param_t& param, params_t... params)
    {
        return sum(func(param),reduce_over_params(func, params...));
    }
}