#pragma once

#include <concepts>
#include <vector>

#include "core/bounding_box.h"
#include "core/ctrs.h"

namespace spade::grid
{
    //This is for GPU kernels that need to efficiently value-capture grid properties
    template <coords::coordinate_system coord_t, const int dim, typename bbox_container_t, typename ibdy_container_t>
    struct grid_geometry_t
    {
        using float_t = typename coord_t::coord_type;
        using iarr_t = ctrs::array<int, 3>;
        
        //Note that this constitutes a local-only view of the geometry of a grid
        coord_t             coords;         // Note that this currently constitutes only an analytical coordinate system
        ctrs::array<int, 3> num_cell;
        ctrs::array<int, 3> num_exch;
        bbox_container_t    bounding_boxes; // Bounding box for each block, comp. coordinates
        ibdy_container_t    domain_boundary;
        
        grid_geometry_t(){}
        template <ctrs::basic_array carr_t>
        grid_geometry_t(const coord_t& coords_in, const carr_t& num_cell_in, const carr_t& num_exch_in)
        : coords{coords_in}
        {
            ctrs::copy_array(num_cell_in, num_cell, 1);
            ctrs::copy_array(num_exch_in, num_exch, 0);
            
            //Geometry data is currently populated outside, need to make this a little cleaner.
        }
        
        const auto& is_domain_boundary(const std::size_t& lb) const { return domain_boundary[lb]; }
        const auto& get_bounding_box  (const std::size_t& lb) const { return bounding_boxes[lb]; }
        const auto& get_coords() const { return coords; }
        const auto& image() const { return *this; }
    };
}
