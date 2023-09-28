#pragma once

#include <concepts>
#include "core/coord_system.h"
#include "grid/unstructured_geometry.h"
namespace spade::grid
{
    //todo: make a concept for the element type
    template
    <
        typename uelem_t,
        coords::coordinate_system coord_t,
        parallel::parallel_group par_group_t
    >
    struct unstructured_grid_t
    {
        using element_type        = uelem_t;
        using dtype               = coord_t::coord_type;
        using coord_type          = coord_t::coord_type;
        using coord_sys_type      = coord_t;
        using coord_point_type    = coords::point_t<coord_type>;
        using group_type          = par_group_t;
        using geometry_type       = unstructured_geometry_t<uelem_t, coord_t, device::shared_vector>;
        using geometry_image_type = unstructured_geometry_t<uelem_t, coord_t, utils::const_raw_ptr_t>; //Going to need something more sophisticated if we want to do FSI
        
        const group_type& grid_group;
        // geometry_type     local_geometry;
        // geometry_type     global_geometry;
        
        //                                                  |
        //                                                  V this should go into the geometry!
        unstructured_grid_t(const uelem_t&, const coord_t& coord_sys_in, const par_group_t& group_in) : grid_group{group_in}
        {
            
        }
    };
}