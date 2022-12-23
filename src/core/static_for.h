#pragma once

#include <concepts>

namespace spade::algs
{
    template <const int imin, const int imax, typename callable_t>
    requires (imin >= imax)
    void static_for(const callable_t& func){}
    
    template <const int imin, const int imax, typename callable_t>
    requires (imin < imax)
    void static_for(const callable_t& func)
    {
        func(std::integral_constant<int, imin>());
        static_for<imin+1, imax>(func);
    }
}