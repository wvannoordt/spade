#pragma once

#include "core/ctrs.h"

namespace spade::utils
{
    template <typename arr_t>
    constexpr static arr_t permute(const arr_t& arr, const ctrs::array<int, arr_t::size()>& perm)
    {
        arr_t output;
        for (std::size_t i = 0; i < arr_t::size(); ++i)
        {
            output[i] = arr[perm[i]];
        }
        return output;
    }
    
    template <typename arr_t>
    constexpr static arr_t ipermute(const arr_t& arr, const ctrs::array<int, arr_t::size()>& perm)
    {
        arr_t output;
        for (std::size_t i = 0; i < arr_t::size(); ++i)
        {
            output[perm[i]] = arr[i];
        }
        return output;
    }
}