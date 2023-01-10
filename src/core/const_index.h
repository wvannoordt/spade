#pragma once
#include <type_traits>
#include "static_math.h"
namespace spade::udci
{
    using integral_t = std::size_t;
    template <integral_t ii>
    struct idx_const_t {constexpr static integral_t value = ii;};

    template <integral_t ii>
    constexpr idx_const_t<ii> idx_const;

    constexpr integral_t c_to_i(char c)
    {
        return c - '0';
    }

    constexpr integral_t constexpr_pow(integral_t base, integral_t exp)
    {
        integral_t ret = 1;
        for (integral_t i = 0; i < exp; ++i) {
            ret *= base;
        }
        return ret;
    }

    template <char... chars_v, integral_t... ints_v>
    constexpr integral_t to_size_t_impl(std::index_sequence<ints_v...>)
    {
        return ((c_to_i(chars_v) * constexpr_pow(10, sizeof...(ints_v) - 1 - ints_v)) + ...);
    }

    template <char... chars_v>
    constexpr integral_t to_size_t()
    {
        return to_size_t_impl<chars_v...>(std::make_index_sequence<sizeof...(chars_v)>{});
    }
}

template <char... chars_v>
constexpr auto operator""_c()
{
    using namespace spade::udci;
    return idx_const<to_size_t<chars_v...>()>;
}