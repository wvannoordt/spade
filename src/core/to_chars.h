#pragma once

#include "core/ctrs.h"
#include "core/unsafe_cast.h"

namespace spade::utils
{
    template <typename data_t>
    requires(std::is_trivially_copyable<data_t>::value)
    _sp_hybrid inline auto to_chars(const data_t& data)
    {
        using output_t = ctrs::array<char, sizeof(data_t)>;
        return algs::unsafe_cast<output_t>(data);
    }
    
    template <typename data_t>
    requires(std::is_trivially_copyable<data_t>::value)
    _sp_hybrid inline data_t from_chars(const ctrs::array<char, sizeof(data_t)>& data)
    {
        return algs::unsafe_cast<data_t>(data);
    }
}