#pragma once

#include "grid/grid.h"

namespace spade::ibm
{
    template <typename grid_t, typename surf_geom_t>
    static void trace_xrays(const grid_t& grid, const surf_geom_t& geom)
    {
        const auto  is_intersect = [&](const auto& lb_loc)
        {
            const auto bnd = grid.get_bounding_box(lb_loc);
            return geom.partially_contained_by(bnd);
        };
        
        const auto  lbs = grid.select_blocks(is_intersect, partition::local);
    }
}