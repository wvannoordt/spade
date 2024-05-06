#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct value : public info_base<value> // interpolation to nodes?
        {
            constexpr static bool requires_direction = false;
            constexpr static bool is_shmem_buffered  = true;

            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename array_t::alias_type;
            
            std::string get_name() const
            {
                return "value";
            }
            
            template <typename grid_view_t, typename array_t, typename index_t>
            requires(index_t::centering_type() == array_t::centering_type())
            _sp_hybrid static void compute(
                const grid_view_t&,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                out = array.get_elem(idx);
            }

            template <typename grid_view_t, typename array_t, typename index_t>
            requires((index_t::centering_type() == grid::face_centered) && (grid::cell_centered == array_t::centering_type()))
            _sp_hybrid static void compute(
                const grid_view_t&,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                using float_t = typename array_t::value_type;
                out =  array.get_elem(grid::face_to_cell(idx,0));
                out += array.get_elem(grid::face_to_cell(idx,1));
                out *= float_t(0.5);
            }
        };
    }
}