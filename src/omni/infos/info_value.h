#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct value : public info_base<value> // interpolation to nodes?
        {
            constexpr static bool requires_direction = false;

            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename array_t::alias_type;
            
            template <typename grid_view_t, typename array_t, typename index_t>
            requires(index_t::centering_type() == array_t::centering_type())
            static void compute(
                const grid_view_t&,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out = array.get_elem(idx);
            }

            template <typename grid_view_t, typename array_t, typename index_t>
            requires((index_t::centering_type() == grid::face_centered) && (grid::cell_centered == array_t::centering_type()))
            static void compute(
                const grid_view_t&,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out =  array.get_elem(grid::face_to_cell(idx,0));
                out += array.get_elem(grid::face_to_cell(idx,1));
                out *= 0.5;
            }
        };
    }
}