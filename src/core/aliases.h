#pragma once

#include <tuple>
#include <vector>

namespace spade::aliases
{
    template <typename... types_t> using vector = std::vector<types_t...>;
    template <typename... types_t> using tuple  = std::tuple<types_t...>;
}