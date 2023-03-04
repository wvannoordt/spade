#pragma once
#include <type_traits>
#include "static_math.h"
namespace spade::udci
{
    using integral_t = int;
    template <integral_t ii>
    struct idx_const_t
    {
        constexpr static integral_t value = ii;
        auto operator -() const {return idx_const_t<-ii>();}
    };
    namespace detail
    {
        constexpr integral_t c_to_i(char c) { return c - '0'; }
        
        template <char v0, char... chars_v> struct compute_b10_val_t
        {
            constexpr static int value = c_to_i(v0)*static_math::pow<10, sizeof...(chars_v)>::value + compute_b10_val_t<chars_v...>::value;
        };
        template <char v0> struct compute_b10_val_t<v0>
        {
            constexpr static int value = c_to_i(v0);
        };
    }
}

template <char... chars_v>
constexpr auto operator""_c()
{
    const int ii = spade::udci::detail::compute_b10_val_t<chars_v...>::value;
    return spade::udci::idx_const_t<ii>();
}