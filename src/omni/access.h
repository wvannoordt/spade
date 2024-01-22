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