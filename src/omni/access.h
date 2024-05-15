#pragma once

#include <concepts>

#include "omni/info.h"
#include "omni/data.h"

#include "omni/buffered.h"

namespace spade::omni
{
    template <typename info_t, typename data_t>
    requires(data_t::list_type::template contains<info_t>)
    _sp_inline
    const typename info_t::array_data_type<typename data_t::array_type, data_t::info_center>&
    _sp_hybrid access(const data_t& data)
    {
        static_assert(data_t::list_type::template contains<info_t>, "attempted static access on info type not contained in stencil data");
        if constexpr (std::same_as<info_t, typename data_t::info_type>)
        {
            return data.data;
        }
        else
        {
            return access<info_t>(data.next);
        }
    }
    
    template <typename info_t, typename data_t>
    requires(data_t::list_type::template contains<info_t>)
    _sp_inline
    typename info_t::array_data_type<typename data_t::array_type, data_t::info_center>&
    _sp_hybrid access(data_t& data)
    {
        // static_assert(data_t::list_type::template contains<info_t>, "attempted static access on info type not contained in stencil data");
        if constexpr (std::same_as<info_t, typename data_t::info_type>)
        {
            return data.data;
        }
        else
        {
            return access<info_t>(data.next);
        }
    }
    
    template <typename info_t, typename offst_t, typename data_t>
    _sp_inline auto& _sp_hybrid access_at(data_t& data)
    {
        constexpr grid::array_centering ctr = offst_t::template relative_node_centering<data_t::stencil_type::center()>;
        constexpr int idx = index_of<typename data_t::stencil_type, offst_t>;
        using idx_t = udci::idx_const_t<idx>;
        return access<info_t>(data.template seek_element<ctr>(idx_t()));
    }
    
    template <typename info_t, typename offst_t, typename data_t>
    _sp_inline const auto& _sp_hybrid access_at(const data_t& data)
    {
        constexpr grid::array_centering ctr = offst_t::template relative_node_centering<data_t::stencil_type::center()>;
        constexpr int idx = index_of<typename data_t::stencil_type, offst_t>;
        using idx_t = udci::idx_const_t<idx>;
        return access<info_t>(data.template seek_element<ctr>(idx_t()));
    }
    
    template <typename info_t, typename local_data_t, typename shared_data_t>
    requires(local_data_t::list_type::template contains<info_t>)
    _sp_inline _sp_hybrid const auto& access(const disparate_data_t<local_data_t, shared_data_t>& data)
    {
        return access<info_t>(data.local);
    }
    
    template <typename info_t, typename local_data_t, typename shared_data_t>
    requires(shared_data_t::list_type::template contains<info_t>)
    _sp_inline _sp_hybrid const auto& access(const disparate_data_t<local_data_t, shared_data_t>& data)
    {
        return access<info_t>(data.shared);
    }
}