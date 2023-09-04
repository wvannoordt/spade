#pragma once

#include "omni/info.h"

namespace spade::omni
{
    template <typename info_t, typename grid_view_t, typename array_t, typename index_t, typename direction_t, typename data_t>
    requires (info_t::requires_direction)
    _sp_hybrid void invoke_compute(const grid_view_t& grid, const array_t& array, const index_t& index, const direction_t& direction, data_t& data)
    {
        info_t::compute(grid, array, index, direction, data);
    }

    template <typename info_t, typename grid_view_t, typename array_t, typename index_t, typename direction_t, typename data_t>
    _sp_hybrid void invoke_compute(const grid_view_t& grid, const array_t& array, const index_t& index, const direction_t& direction, data_t& data)
    {
        info_t::compute(grid, array, index, data);
    }

    template <typename grid_view_t, typename index_t, typename array_t, typename direction_t, typename stencil_data_t>
    _sp_hybrid void retrieve_impl(const grid_view_t& grid, const array_t& array, const index_t& idx_center, const direction_t& direction, stencil_data_t& data)
    {
        const int num_elem = stencil_data_t::stencil_type::num_elements();
        algs::static_for<0,num_elem>([&](const auto& ii)
        {
            const int i = ii.value;
            
            //fetch the list of info that we need at the current stencil element
            auto& info_list = detail::get_info_list_at<i,0>(data);
            using offset_type = typename stencil_data_t::stencil_type::stencil_element<i>::position_type;
            const auto idx  = offset_type::compute_index(idx_center);
            algs::static_for<0,info_list.num_info_elems()>([&](const auto& jj)
            {
                const int j = jj.value;
                auto& info_list_data = detail::get_info_list_data_at<j,0>(info_list);
                using info_type = typename utils::remove_all<decltype(info_list_data)>::type::info_type;
                // constexpr bool use_dir = requires{info_type::compute(grid, array, idx, direction, info_list_data.data);};
                invoke_compute<info_type>(grid, array, idx, direction, info_list_data.data);
            });
        });
    }

    template <typename grid_view_t, typename index_t, typename array_t, typename stencil_data_t>
    requires(requires{index_t().dir();})
    _sp_hybrid void retrieve(const grid_view_t& grid, const array_t& array, const index_t& idx_center, stencil_data_t& data)
    {
        retrieve_impl(grid, array, idx_center, idx_center.dir(), data);
    }

    template <typename grid_view_t, typename index_t, typename array_t, typename stencil_data_t>
    _sp_hybrid void retrieve(const grid_view_t& grid, const array_t& array, const index_t& idx_center, stencil_data_t& data)
    {
        retrieve_impl(grid, array, idx_center, info::undirected_tag_t(), data);
    }
}