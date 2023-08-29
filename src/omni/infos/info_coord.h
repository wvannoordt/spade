#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct coord : public info_base<coord>
        {
            constexpr static bool requires_direction = false;
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename array_t::grid_type::coord_point_type;
            
            template <typename grid_view_t, typename array_t, typename index_t>
            _sp_hybrid static void compute(
                const grid_view_t& grid,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out = grid.get_coords(idx);
            }
        };
    }
}