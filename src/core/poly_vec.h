#pragma once

#include <vector>

#include "core/utils.h"

namespace spade::utils
{
    template <typename data_t, typename... datas_t>
    requires(is_unique<data_t, datas_t...>)
    struct poly_vec_t
    {
        device::shared_vector<data_t> data;
        poly_vec_t<datas_t...> sub;
        
        template <typename new_t>
        void add(const new_t& item)
        {
            if constexpr (std::same_as<new_t, data_t>)
            {
                data.push_back(item);
            }
            else
            {
                sub.add(item);
            }
        }
        
        template <typename op_t>
        void foreach(const op_t& op) const
        {
            for (const auto& ii: data) op(ii);
            sub.foreach(op);
        }
        
        template <const int ii>
        requires (ii >= 0)
        const auto& operator[] (const udci::idx_const_t<ii>&) const
        {
            if constexpr (ii == 0) return data;
            else return sub[udci::idx_const_t<ii-1>()];
        }
        
        template <const int ii>
        requires (ii >= 0)
        auto& operator[] (const udci::idx_const_t<ii>&)
        {
            if constexpr (ii == 0) return data;
            else return sub[udci::idx_const_t<ii-1>()];
        }
    };
    
    template <typename data_t> struct poly_vec_t<data_t>
    {
        device::shared_vector<data_t> data;
        void add(const data_t& item)
        {
            data.push_back(item);
        }
        
        template <typename op_t>
        void foreach(const op_t& op) const
        {
            for (const auto& ii: data) op(ii);
        }
        
        template <const int ii>
        requires (ii == 0)
        const auto& operator[] (const udci::idx_const_t<ii>&) const
        {
            return data;
        }
        
        template <const int ii>
        requires (ii == 0)
        auto& operator[] (const udci::idx_const_t<ii>&)
        {
            return data;
        }
    };
}