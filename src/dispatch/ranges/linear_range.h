#pragma once

#include "core/ctrs.h"

namespace spade::dispatch::ranges
{
    template <typename idx_t, typename r_device_t>
    struct linear_range_t
    {
        using index_type = ctrs::array<idx_t, 1>;
        using device_t   = r_device_t;
        
        r_device_t dvc;
        idx_t lower;
        idx_t upper;
        
        linear_range_t(const idx_t& lw, const idx_t& up, const r_device_t& de) : dvc{de}, lower{lw}, upper{up} {}
        
        device_t device() const { return dvc; }
        
        bound_box_t<idx_t, 1> idx_bound_box() const
        {
            bound_box_t<idx_t, 1> output;
            output.min(0) = lower;
            output.max(0) = upper;
            return output;
        }
    };
    
    template <ctrs::basic_array arr_t, typename device_t>
    static auto from_array(const arr_t& arr, const device_t& device)
    {
        return linear_range_t(0UL*arr.size(), arr.size(), device);
    }
}