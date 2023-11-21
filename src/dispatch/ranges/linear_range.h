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
    
    // to-do: get rid of this crap
    template <typename idx_t, typename r_device_t>
    struct bilinear_range_t
    {
        using index_type = ctrs::array<idx_t, 2>;
        using device_t   = r_device_t;
        
        r_device_t dvc;
        idx_t i0, i1;
        idx_t j0, j1;
        
        bilinear_range_t(const idx_t& i0_in, const idx_t& i1_in, const idx_t& j0_in, const idx_t& j1_in, const r_device_t& de) : dvc{de}, i0{i0_in}, i1{i1_in}, j0{j0_in}, j1{j1_in} {}
        
        device_t device() const { return dvc; }
        
        bound_box_t<idx_t, 2> idx_bound_box() const
        {
            bound_box_t<idx_t, 2> output;
            output.min(0) = i0;
            output.max(0) = i1;
            output.min(1) = j0;
            output.max(1) = j1;
            return output;
        }
    };
}