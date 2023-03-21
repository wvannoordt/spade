#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct normal : public info_base<normal> // only supported at locations centered at faces, expressed in index space
        {
            constexpr static bool supports_undirected = false;
            constexpr static bool supports_directed   = true;
            template <typename array_t, const grid::array_centering center>
            using array_data_type = ctrs::array<int, array_t::grid_type::coord_point_type::size()>;
            
            template <typename array_t, typename index_t>
            static void compute(
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                print("normal");
            }
        };
    }
}