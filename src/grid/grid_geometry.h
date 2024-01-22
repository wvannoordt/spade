#pragma once

#include <concepts>
#include <vector>

#include "core/bounding_box.h"
#include "core/ctrs.h"

namespace spade::grid
{
    //This is for GPU kernels that need to efficiently value-capture grid properties
    template <
        coords::coordinate_system coord_t,
        const int gdim,
        template <typename> typename container_tt>
    struct grid_geometry_t
    {
        using coord_type = typename coord_t::coord_type;
        using float_t    = typename coord_t::coord_type;
        using point_type = coords::point_t<float_t>;
        using iarr_t     = ctrs::array<int, 3>;
        
        //Note that this constitutes a local-only view of the geometry of a grid
        coord_t                               coords;         // Note that this currently constitutes only an analytical coordinate system
        ctrs::array<int, 3>                   num_cell;
        container_tt<bound_box_t<float_t, 3>> bounding_boxes; // Bounding box for each block, comp. coordinates
        container_tt<bound_box_t<bool, 3>>    domain_boundary;
        ctrs::array<container_tt<int>, 6>     boundary_blocks;
        
        constexpr static int dim() { return gdim; }
        
        grid_geometry_t(){}
        
        template <ctrs::basic_array carr_t>
        grid_geometry_t(const coord_t& coords_in, const carr_t& num_cell_in)
        : coords{coords_in}
        {
            ctrs::copy_array(num_cell_in, num_cell, 1);
            
            //Geometry data is currently populated outside, need to make this a little cleaner.
        }
        
        _sp_hybrid const auto& is_domain_boundary(const std::size_t& lb) const { return domain_boundary[lb]; }
        _sp_hybrid const auto& get_bounding_box  (const std::size_t& lb) const { return bounding_boxes[lb];  }
        _sp_hybrid const auto& get_coords()    const { return coords; }
        _sp_hybrid const auto& get_coord_sys() const { return coords; }
        
        template <typename idx_t>
        _finline_ _sp_hybrid point_type get_coords(const idx_t& i) const
        {
            return coords.map(this->get_comp_coords(i));
        }
        
        template <typename idx_t>
        _finline_ _sp_hybrid point_type get_comp_coords(const idx_t& i) const
        {
            const auto& bbx   = this->get_bounding_box(i.lb());
            const auto  idx_r = get_index_coord(i);
            point_type output;
            if constexpr (dim() == 2) output[2] = 0.0;
            const auto dx = this->get_dx(i.lb());
            for (int d = 0; d < dim(); ++d)
            {
                output[d] = bbx.min(d)+idx_r[d]*dx[d];
            }
            return output;
        }
        
        _sp_hybrid ctrs::array<float_t, 3> get_dx(const std::size_t& lb) const
        {
            const auto& bbx = this->get_bounding_box(lb);
            ctrs::array<float_t, 3> output;
            if constexpr (dim() == 2) output[2] = 1.0;
            for (int d = 0; d < dim(); ++d)
            {
                output[d] = bbx.size(d)/num_cell[d];
            }
            return output;
        }
        
        _sp_hybrid float_t get_dx(const int i, const std::size_t& lb)  const
        {
            return this->get_bounding_box(lb).size(i)/num_cell[i];
        }
    };
}
