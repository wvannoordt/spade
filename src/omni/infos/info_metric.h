#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct metric : public info_base<metric>
        {
            constexpr static bool requires_direction = true;

            template <typename array_t, const grid::array_centering center>
            using array_data_type = ctrs::array<typename array_t::grid_type::coord_type,3>;
            
            template <typename array_t, typename index_t>
            static void compute(
                const array_t& array,
                const index_t& idx,
                const int& dir,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                const auto& grid = array.get_grid();
                out = calc_normal_vector(grid.get_coord_sys(), grid.get_coords(idx), idx, dir);
            }
        };
    }
}