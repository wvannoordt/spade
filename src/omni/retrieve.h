#pragma once

#include "omni/info.h"

namespace spade::omni
{
    template <typename info_t, typename grid_view_t, typename array_t, typename index_t, typename direction_t, typename data_t>
    requires (info_t::requires_direction)
    _sp_inline _sp_hybrid void invoke_compute(const grid_view_t& grid, const array_t& array, const index_t& index, const direction_t& direction, data_t& data)
    {
        info_t::compute(grid, array, index, direction, data);
    }

    template <typename info_t, typename grid_view_t, typename array_t, typename index_t, typename direction_t, typename data_t>
    _sp_inline _sp_hybrid void invoke_compute(const grid_view_t& grid, const array_t& array, const index_t& index, const direction_t& direction, data_t& data)
    {
        info_t::compute(grid, array, index, data);
    }

    template <typename grid_view_t, typename index_t, typename array_t, typename direction_t, typename stencil_data_t>
    _sp_inline _sp_hybrid void retrieve_impl(const grid_view_t& grid, const array_t& array, const index_t& idx_center, const direction_t& direction, stencil_data_t& data)
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
    _sp_inline _sp_hybrid void retrieve(const grid_view_t& grid, const array_t& array, const index_t& idx_center, stencil_data_t& data)
    {
        retrieve_impl(grid, array, idx_center, idx_center.dir(), data);
    }

    template <typename grid_view_t, typename index_t, typename array_t, typename stencil_data_t>
    _sp_inline _sp_hybrid void retrieve(const grid_view_t& grid, const array_t& array, const index_t& idx_center, stencil_data_t& data)
    {
        retrieve_impl(grid, array, idx_center, info::undirected_tag_t(), data);
    }
    
    //Shared memory, buffered version
    template <typename grid_view_t, typename index_t, typename array_t, typename stencil_data_t, typename shmem_t>
    _sp_inline _sp_hybrid void retrieve_shared(const grid_view_t& grid, const array_t& array, const index_t& idx_center, stencil_data_t& data, shmem_t& buffer)
    {
        using stencil_type = typename stencil_data_t::stencil_type;
        constexpr static int face_idx = index_of<stencil_type, offset_t< 0, 0, 0>>;
        constexpr static int cell_idx = index_of<stencil_type, offset_t< 1, 0, 0>>;
        const int face_buf_idx = data.get_face_buf_index(udci::idx_const_t<face_idx>());
        const int cell_buf_idx = data.get_cell_buf_index(udci::idx_const_t<cell_idx>());
        
        using i_c = typename stencil_data_t::cell_data_location;
        using i_f = typename stencil_data_t::face_data_location;
        
        auto& c_data = buffer[i_c()][cell_buf_idx];
        auto& f_data = buffer[i_f()][face_buf_idx];
        
        using c_info_list_type = stencil_data_t::cell_data_type;
        using f_info_list_type = stencil_data_t::face_data_type;
        
        //Buffer cell info
        auto cell_glob_idx = offset_t<1, 0, 0>::compute_index(idx_center);
        algs::static_for<0, c_info_list_type::num_info_elems()>([&](const auto& jj)
        {
            const int j = jj.value;
            auto& info_list_data = detail::get_info_list_data_at<j,0>(c_data);
            using info_type = typename utils::remove_all<decltype(info_list_data)>::type::info_type;
            invoke_compute<info_type>(grid, array, cell_glob_idx, idx_center.dir(), info_list_data.data);
        });
        
        constexpr int num_suppl_l = stencil_data_t::n_l_cells;
        constexpr int num_suppl_r = stencil_data_t::n_r_cells - 1;
        
        int nface = buffer[i_f()].size();
        int num_extra = num_suppl_l + num_suppl_r;
        if (data.line_idx < num_extra)
        {
            int new_buf_idx = cell_buf_idx - num_suppl_l;
            if (data.line_idx >= num_suppl_l) new_buf_idx = data.line_idx + nface;
            auto& c_data_extra = buffer[i_c()][new_buf_idx];
            cell_glob_idx.i(idx_center.dir()) = new_buf_idx - num_suppl_l;
            algs::static_for<0, c_info_list_type::num_info_elems()>([&](const auto& jj)
            {
                const int j = jj.value;
                auto& info_list_data = detail::get_info_list_data_at<j,0>(c_data_extra);
                using info_type = typename utils::remove_all<decltype(info_list_data)>::type::info_type;
                invoke_compute<info_type>(grid, array, cell_glob_idx, idx_center.dir(), info_list_data.data);
            });
        }
        
        //Buffer face info
        algs::static_for<0, f_info_list_type::num_info_elems()>([&](const auto& jj)
        {
            const int j = jj.value;
            auto& info_list_data = detail::get_info_list_data_at<j,0>(f_data);
            using info_type = typename utils::remove_all<decltype(info_list_data)>::type::info_type;
            invoke_compute<info_type>(grid, array, idx_center, idx_center.dir(), info_list_data.data);
        });
    }
}