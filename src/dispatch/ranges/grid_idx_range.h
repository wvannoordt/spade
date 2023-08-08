#pragma once

#include "grid/grid_index_types.h"

namespace spade::dispatch::ranges
{
    template <grid::grid_index index_t, typename r_device_t>
    struct grid_idx_range_t
    {
        using index_type = index_t;
        using device_t   = r_device_t;
        index_t lower;
        index_t upper;
        grid_idx_range_t(const index_t& ll, const index_t& ur, const r_device_t& device)
        : lower{ll}, upper{ur} {}
        
        template <typename kernel_t>
        void execute(const kernel_t kernel)
        {
            //do stuff
        }
    };
}