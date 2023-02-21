#pragma once

#include "core/grid_array.h"

#include "omni/geom.h"

namespace spade::omni
{
    struct end_data_t {};
    
    template <const int i0, const int i1, typename info_list_t, typename array_t, const grid::array_centering center>
    struct info_list_data_t
    {
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
        
        template <udci::integral_t ii> const auto& cell(const udci::idx_const_t<ii>& idx) {return next;};
        template <udci::integral_t ii> const auto& face(const udci::idx_const_t<ii>& idx) {return next;};
        template <udci::integral_t ii> const auto& node(const udci::idx_const_t<ii>& idx) {return next;};
        template <udci::integral_t ii> const auto& edge(const udci::idx_const_t<ii>& idx) {return next;};
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