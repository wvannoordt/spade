#pragma once

#include "core/bounding_box.h"
#include "grid/grid_index_types.h"

namespace spade::dispatch::ranges
{
    template <grid::grid_index index_t, typename r_device_t>
    struct grid_idx_range_t
    {
        using index_type = index_t;
        using device_t   = r_device_t;
        r_device_t dvc;
        index_t lower;
        index_t upper;
        grid_idx_range_t(const index_t& ll, const index_t& ur, const r_device_t& device)
        : lower{ll}, upper{ur}, dvc{device} {}
        
        device_t device() const { return dvc; }
        
        bound_box_t<int, index_t::size()> idx_bound_box() const
        {
            bound_box_t<int, index_t::size()> output;
            for (int i = 0; i < index_t::size(); ++i)
            {
                output.min(i) = lower[i];
                output.max(i) = upper[i];
            }
            return output;
        }
    };
}