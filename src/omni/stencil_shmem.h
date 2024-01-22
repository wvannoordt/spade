#pragma once

#include "dispatch/shmem.h"
#include "omni/stencil_union.h"

namespace spade::omni
{
    template <typename omni_type, typename array_t>
    requires (omni_type::center() == grid::face_centered)
    inline auto get_shmem(const std::size_t& num_faces, const array_t&)
    {
        // THIS WILL NEED TO BE FIXED!!!!!
        // This will cause someone a lot of grief at some point. I actually need to combine all cell infos into a single list.
        using full_face_info_list = typename omni_type::info_at<omni::offset_t<0, 0, 0>>;
        using full_cell_info_list = typename omni_type::info_at<omni::offset_t<1, 0, 0>>;
        
        using face_info_list = typename omni::shmem_info_list<full_face_info_list>;
        using cell_info_list = typename omni::shmem_info_list<full_cell_info_list>;
        
        using face_data_type = typename omni::get_list_data<face_info_list, array_t, grid::face_centered>;
        using cell_data_type = typename omni::get_list_data<cell_info_list, array_t, grid::cell_centered>;
        
        const std::size_t left_extent  = -static_math::moddiv<omni_type::template min_extent<0>, 2>::value;
        const std::size_t right_extent = static_math::moddiv<omni_type::template max_extent<0>, 2>::value + 1;
        const std::size_t num_cells   = num_faces - 1 + left_extent + right_extent;
        
        return dispatch::shmem::make_shmem(dispatch::shmem::vec<face_data_type>(num_faces), dispatch::shmem::vec<cell_data_type>(num_cells));
    }
}