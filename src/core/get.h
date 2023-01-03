#pragma once

#include "core/aliases.h"

namespace spade::ctrs
{
    template <const int i, typename tuple_t> const auto& get(const tuple_t& tup)
    {
        if constexpr (requires {tup.get_tuple();}) return std::get<i>(tup.get_tuple());
        else return std::get<i>(tup);
    }
}