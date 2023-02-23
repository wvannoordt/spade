#pragma once

#include "core/grid_array.h"

#include "omni/geom.h"

namespace spade::omni
{
    struct end_data_t {};
    
    template <const int i0, const int i1, typename info_list_t, typename array_t, const grid::array_centering center>
    struct info_list_data_t
    {
        constexpr static grid::array_centering info_center = center;
        using array_type = array_t;
        using info_type = typename info_list_t::info_elem<i0>;
        using data_type = typename info_type::array_data_type<array_t, center>;
        data_type data;
        info_list_data_t<i0+1, i1, info_list_t, array_t, center> next;
    };
    
    template <const int i0, typename info_list_t, typename array_t, const grid::array_centering center>
    struct info_list_data_t<i0, i0, info_list_t, array_t, center>
    {
        
    };
    
    template <const int i0, const int i1, typename stencil_t, typename array_t>
    struct stencil_data_sublist_t
    {
        using stencil_type   = stencil_t;
        using array_type     = array_t;
        using stencil_info_t = typename stencil_t::stencil_element<i0>::info_type;
        using stencil_ofst_t = typename stencil_t::stencil_element<i0>::position_type;
        using next_type      = stencil_data_sublist_t<i0+1, i1, stencil_t, array_t>;
        
        constexpr static grid::array_centering centering_at_i0_elem = stencil_ofst_t::template relative_node_centering<stencil_t::center()>;
        
        info_list_data_t<0, stencil_info_t::num_infos(), stencil_info_t, array_t, centering_at_i0_elem> data;
        next_type next;
        
        //const qualified
        template <const grid::array_centering ctr, udci::integral_t ii>
        const auto& seek_element(const udci::idx_const_t<ii>& idx) const
        {
            if constexpr(ii == 0 && ctr == centering_at_i0_elem)
            {
                return data;
            }
            else if constexpr (ctr == centering_at_i0_elem)
            {
                return next.template seek_element<ctr>(udci::idx_const_t<ii-1>());
            }
            else
            {
                return next.template seek_element<ctr>(udci::idx_const_t<ii>());
            }
        }
        
        //not const qualified
        template <const grid::array_centering ctr, udci::integral_t ii>
        auto& seek_element(const udci::idx_const_t<ii>& idx)
        {
            if constexpr(ii == 0 && ctr == centering_at_i0_elem)
            {
                return data;
            }
            else if constexpr (ctr == centering_at_i0_elem)
            {
                return next.template seek_element<ctr>(udci::idx_const_t<ii-1>());
            }
            else
            {
                return next.template seek_element<ctr>(udci::idx_const_t<ii>());
            }
        }
        
        
        //const qualified
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_cell())
        constexpr const auto& cell(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::cell_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_face())
        constexpr const auto& face(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::face_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_node())
        constexpr const auto& node(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::node_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_edge())
        constexpr const auto& edge(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::edge_centered>(idx);
        }
        
        //not const qualified
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_cell())
        constexpr auto& cell(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::cell_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_face())
        constexpr auto& face(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::face_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_node())
        constexpr auto& node(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::node_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < stencil_t::num_edge())
        constexpr auto& edge(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::edge_centered>(idx);
        }
    };
    
    template <const int i0, typename stencil_t, typename array_t>
    struct stencil_data_sublist_t<i0, i0, stencil_t, array_t>
    {
        using stencil_type = stencil_t;
        using array_type   = array_t;
    };

    template <typename stencil_t, typename array_t>
    using stencil_data_t = stencil_data_sublist_t<0, stencil_t::num_elements(), stencil_t, array_t>;
}