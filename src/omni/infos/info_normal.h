#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct normal : public info_base<normal> // only supported at locations centered at faces, expressed in index space
        {
            constexpr static bool requires_direction = true;
            
            template <typename array_t, const grid::array_centering center>
            using array_data_type = ctrs::array<int, array_t::grid_type::coord_point_type::size()>;
            
            template <typename grid_view_t, typename array_t, typename index_t>
            _sp_hybrid static void compute(
                const grid_view_t&,
                const array_t& array,
                const index_t& idx,
                const int& dir,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out = 0;
                out[dir] = 1;
            }
        };
    }
}