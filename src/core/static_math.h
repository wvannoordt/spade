#pragma once

#include <concepts>
#include <type_traits>

namespace spade::static_math
{
    template <const int i1> struct int_const_t { constexpr static int value = i1; };
    
    template <std::floating_point data_t> struct real_const_t
    {
        data_t data;
        real_const_t(const data_t& data_in) {data = data_in;}
    };
    
    template <const int i1, const unsigned int i2> struct pow
    {
        constexpr static int value = 
            std::conditional<(i2==0), int_const_t<1>, int_const_t<i1>>::type::value*
            std::conditional<(i2<=1), int_const_t<1>, pow<i1, i2-1>>::type::value;
    };
}