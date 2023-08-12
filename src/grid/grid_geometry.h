#pragma once

#include <concepts>

#include "core/bounding_box.h"
#include "core/ctrs.h"
/*
namespace spade::grid
{
    template <coords::coordinate_system coord_t, const int dim>
    struct grid_geometry_t
    {
        using float_t = typename coord_t::coord_type;
        
        //Note that this constitutes a local-only view of the geometry of a grid
        coord_t                              coord_system;   // Note that this currently constitutes only an analytical coordinate system
        something_t<bound_box_t<float_t, 3>> bounding_boxes; // Bounding box for each block
        
    };
}
*/