#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct index : public info_base<index>
        {
            constexpr static bool requires_direction = false;
            constexpr static bool is_shmem_buffered  = false;
            
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
            
            template <typename grid_view_t, typename array_t, typename index_t>
            _sp_hybrid static void compute(
                const grid_view_t&,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out = idx;
            }
        };
    }
}