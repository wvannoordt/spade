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
            
            template <typename array_t, typename index_t>
            static void compute(
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out = array.get_grid().get_coords(idx);
            }
        };
    }
}