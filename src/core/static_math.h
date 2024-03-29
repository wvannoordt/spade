#pragma once


#include <type_traits>

#include "core/config.h"

namespace spade::static_math
{
    template <const int i1> struct int_const_t { constexpr static int value = i1; };
    
    template <std::floating_point data_t> struct real_const_t
    {
        data_t data;
        real_const_t(const data_t& data_in) {data = data_in;}
    };
    
    template <const int n> struct factorial
    {
        constexpr static int value = n*factorial<n-1>::value;
    };
    template <> struct factorial<0>
    {
        constexpr static int value = 1;
    };
    
    template <const int i1, const unsigned int i2> struct pow
    {
        constexpr static int value = 
            std::conditional<(i2==0), int_const_t<1>, int_const_t<i1>>::type::value*
            std::conditional<(i2<=1), int_const_t<1>, pow<i1, i2-1>>::type::value;
    };
    
    
    //Note that mod and moddiv must be used together to maintain the invariant that
    // M = K*(M/K) + M mod K
    template <const int val, const int modu> struct mod
    {
        constexpr static int baseval = (val - modu*(val/modu));
        constexpr static int value = (baseval<0)?(baseval+modu):baseval;
    };
    template <const int val, const int div> struct moddiv
    {
        constexpr static int value = (val - ((val<0)?(div-1):0))/div;
    };
}