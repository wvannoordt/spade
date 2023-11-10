#pragma once

#include "sym/symbol.h"

namespace spade::sym::literals
{
    template <typename T, T... vals>
    constexpr auto operator""_sym() { return symbol_t<T, vals...>(); }
}