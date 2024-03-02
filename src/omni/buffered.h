#pragma once

#include "core/static_math.h"

#include "omni/stencil.h"
#include "omni/info.h"
#include "omni/stencil_union.h"

namespace spade::omni
{
    namespace detail
    {
        template <typename accum_t, typename... infos_t>
        struct lremove_buffered_impl;
        
        template <typename accum_t, typename info_t, typename... infos_t>
        struct lremove_buffered_impl<accum_t, info_list_t<info_t, infos_t...>>
        {
            using new_accum_t = typename std::conditional<info_t::is_shmem_buffered, accum_t, typename accum_t::expand_list<info_t>>::type;
            using type        = typename lremove_buffered_impl<new_accum_t, info_list_t<infos_t...>>::type;
        };
        
        template <typename accum_t, typename info_t>
        struct lremove_buffered_impl<accum_t, info_list_t<info_t>>
        {
            using new_accum_t = typename std::conditional<info_t::is_shmem_buffered, accum_t, typename accum_t::expand_list<info_t>>::type;
            using type = new_accum_t;
        };
        
        template <typename thing_t>
        struct remove_buffered_impl;
        
        template <typename... infos_t>
        struct remove_buffered_impl<info_list_t<infos_t...>>
        {
            using type = typename lremove_buffered_impl<info_list_t<>, info_list_t<infos_t...>>::type;
        };
        
        template <const grid::array_centering center_val, typename offset_t, typename list_t, typename... elems_t>
        struct remove_buffered_impl<stencil_t<center_val, elem_t<offset_t, list_t>, elems_t...>>
        {
            using here_type = stencil_t<center_val, elem_t<offset_t, typename remove_buffered_impl<list_t>::type>>;
            using othe_type = remove_buffered_impl<stencil_t<center_val, elems_t...>>::type;
            using type = stencil_union<here_type, othe_type>;
        };
        
        template <const grid::array_centering center_val, typename offset_t, typename list_t>
        struct remove_buffered_impl<stencil_t<center_val, elem_t<offset_t, list_t>>>
        {
            using here_type = stencil_t<center_val, elem_t<offset_t, typename remove_buffered_impl<list_t>::type>>;
            using type = here_type;
        };
    }
    
    
    template <typename thing_t> using remove_buffered = typename detail::remove_buffered_impl<thing_t>::type;
    
    
    template <typename local_t, typename shared_t>
    struct disparate_data_t
    {
        const local_t&  local;
        const shared_t& shared;
    };
    
    template <typename omni_t, typename arr_t, typename shmem_t>
    requires(
        omni_t::template max_extent<0> == -1*omni_t::template min_extent<0>
        && omni_t::template max_extent<1> == 0
        && omni_t::template min_extent<1> == 0
        && omni_t::template min_extent<2> == 0
        && omni_t::template max_extent<2> == 0
    )
    struct buffered_t
    {
        shmem_t& shmem;
        int line_idx;
        int face_idx_base;
        int face_pitch;
        int cell_idx_base;
        int cell_pitch;
        
        using empty0_t = stencil_t<grid::face_centered, elem_t<offset_t<0,0,0>, info_list_t<>>>;
        using empty1_t = stencil_t<grid::face_centered, elem_t<offset_t<1,0,0>, info_list_t<>>>;
        
        
        using stencil_type = stencil_union<omni_t, empty0_t, empty1_t>;
        using omni_type    = remove_buffered<omni_t>;
        using data_type    = stencil_data_t<omni_type, arr_t>;
        
        using face_data_location = udci::idx_const_t<0>;
        using cell_data_location = udci::idx_const_t<1>;
        
        using full_face_info_list = typename omni_t::info_at<omni::offset_t<0, 0, 0>>;
        using full_cell_info_list = typename omni_t::info_at<omni::offset_t<1, 0, 0>>;
        
        using face_info_list = typename omni::shmem_info_list<full_face_info_list>;
        using cell_info_list = typename omni::shmem_info_list<full_cell_info_list>;
        
        using face_data_type = typename omni::get_list_data<face_info_list, arr_t, grid::face_centered>;
        using cell_data_type = typename omni::get_list_data<cell_info_list, arr_t, grid::cell_centered>;
        
        constexpr static int n_l_cells  = -static_math::moddiv<stencil_type::template min_extent<0>,2>::value;
        constexpr static int n_r_cells  =  static_math::moddiv<stencil_type::template max_extent<0>,2>::value + 1;
        
        data_type supplemental;
        
        template <const int ii>
        _sp_hybrid int get_cell_buf_index(const udci::idx_const_t<ii>&) const
        {
            using offset_type = offset_at<stencil_type, grid::cell_centered, ii>;
            constexpr static int cell_offst = static_math::moddiv<offset_type::template elem<0>,2>::value;
            int ordinary = line_idx;
            return cell_idx_base + cell_pitch*ordinary;
        }
        
        template <const int ii>
        _sp_hybrid int get_face_buf_index(const udci::idx_const_t<ii>&) const
        {
            using offset_type = offset_at<stencil_type, grid::face_centered, ii>;
            int ordinary = line_idx;
            return face_idx_base + face_pitch*ordinary;
        }
        
        //const qualified
        template <udci::integral_t ii>
        requires (ii < stencil_type::num_cell())
        _sp_hybrid constexpr const auto cell(const udci::idx_const_t<ii>& idx) const
        {
            int buf_index = get_cell_buf_index(idx);
            const auto& shmem_data = shmem[cell_data_location()][cell_idx_base + cell_pitch*buf_index];
            const auto& local_data = supplemental.cell(idx);
            using shared_t = typename utils::remove_all<decltype(shmem_data)>::type;
            using local_t = typename utils::remove_all<decltype(local_data)>::type;
            return disparate_data_t<local_t, shared_t>{local_data, shmem_data};
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_type::num_face())
        _sp_hybrid constexpr const auto face(const udci::idx_const_t<ii>& idx) const
        {
            int buf_index = get_face_buf_index(idx);
            const auto& shmem_data = shmem[face_data_location()][face_idx_base + face_pitch*buf_index];
            const auto& local_data = supplemental.face(idx);
            using shared_t = typename utils::remove_all<decltype(shmem_data)>::type;
            using local_t =  typename utils::remove_all<decltype(local_data)>::type;
            return disparate_data_t<local_t, shared_t>{local_data, shmem_data};
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_type::num_node())
        _sp_hybrid constexpr const auto& node(const udci::idx_const_t<ii>& idx) const
        {
            static_assert(stencil_type::num_node()==0, "Nodes not implemented for omni shmem types");
            return supplemental.node(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_type::num_edge())
        _sp_hybrid constexpr const auto& edge(const udci::idx_const_t<ii>& idx) const
        {
            static_assert(stencil_type::num_edge()==0, "Edges not implemented for omni shmem types");
            return supplemental.edge(idx);
        }
        
        _sp_hybrid const auto root() const
        {
            constexpr int idx_v = index_of<stencil_type, offset_t<0,0,0>>;
            return this->face(udci::idx_const_t<idx_v>());
        }
        
        //const qualified
        template <const grid::array_centering ctr, udci::integral_t ii>
        requires (ctr == grid::face_centered)
        _sp_hybrid const auto seek_element(const udci::idx_const_t<ii>& idx) const
        {
            return this->face(idx);
        }
        
        template <const grid::array_centering ctr, udci::integral_t ii>
        requires (ctr == grid::cell_centered)
        _sp_hybrid const auto seek_element(const udci::idx_const_t<ii>& idx) const
        {
            return this->cell(idx);
        }
    };
    
    template <typename istencil_t, typename array_t, typename shmem_t>
    _sp_hybrid auto make_buffered_data(shmem_t& shmem, int i_face_line, int face_base_idx, int face_pitch, int cell_base_idx, int cell_pitch)
    {
        return buffered_t<istencil_t, array_t, shmem_t>{shmem, i_face_line, face_base_idx, face_pitch, cell_base_idx, cell_pitch};
    }
}