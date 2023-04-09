#pragma once

#include "omni/stencil_union.h"

namespace spade::omni
{
    //Used to create objects with the call operator ()(const auto&)
    //from a bunch of objects that don't 
    template <typename callable_t, typename... omnis_t>
    struct omni_kernel_t
    {
        const callable_t& func;
        const std::tuple<omnis_t...> args;
        using omni_type = combine_omni_stencils<omnis_t...>;
        omni_kernel_t(const callable_t& func_in, const omnis_t&... stuff) : func{func_in}, args(std::make_tuple(stuff...)) {}
        auto operator () (const auto& input) const
        {
            return std::apply([&](const auto&... stuff){return func(stuff..., input);}, args);
        }
    };

    template <typename callable_t, typename... omnis_t>
    static auto make_kernel(const callable_t& func, const omnis_t&... omnis)
    {
        return omni_kernel_t(func, omnis...);
    }
}