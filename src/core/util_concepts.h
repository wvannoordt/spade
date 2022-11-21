#pragma once


#include "core/config.h"
namespace spade::utils
{
    template <typename func_t, typename... args_t> concept const_callable = requires(const func_t& f, const args_t&... args)
    {
        f(args...);
    };
}