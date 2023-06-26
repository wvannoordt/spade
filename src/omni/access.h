#pragma once

#include <concepts>

#include "omni/info.h"
#include "omni/data.h"

namespace spade::omni
{
    template <typename info_t, typename data_t>
    const typename info_t::array_data_type<typename data_t::array_type, data_t::info_center>&
    access(const data_t& data)
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
    typename info_t::array_data_type<typename data_t::array_type, data_t::info_center>&
    access(data_t& data)
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
}