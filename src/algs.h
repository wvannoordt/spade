#pragma once

namespace cvdf::algs
{
    template <class ltype, class rtype> void sum_reduce(ltype& total, const rtype& r)
    {
        total += r;
    }
    
    template <class reduction_t, class invokable_t, class param_t>
    static auto reduce_over_params(const reduction_t& reduce, const invokable_t& func, const param_t& param)
    {
        return func(param);
    }
    
    template <class reduction_t, class invokable_t, class param_t, class... params_t>
    static auto reduce_over_params(const reduction_t& reduce_op, const invokable_t& func, const param_t& param, params_t... params)
    {
        auto total = func(param);
        reduce_op(total, reduce_over_params(reduce_op, func, params...));
        return total;
    }
    template <class invokable_t, class param_t, class... params_t>
    static auto reduce_over_params(const invokable_t& func, const param_t& param, params_t... params)
    {
        typedef decltype(func(param)) ret_type;
        return reduce_over_params(sum_reduce<ret_type,ret_type>, func, param, params...);
    }
    
    template <class invokable_t, class param_t>
    static void foreach_param(const invokable_t& func, const param_t& param)
    {
        func(param);
    }
    
    template <class invokable_t, class param_t, class... params_t>
    static void foreach_param(const invokable_t& func, const param_t& param, params_t... params)
    {
        func(param);
        foreach_param(params...);
    }
    
    
}