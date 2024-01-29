#pragma once

#include <syncstream>

namespace spade::grid
{
    namespace detail
    {
        template <class T> concept has_value_type = requires(T t) {typename T::value_type;};
        template <typename data_t> struct get_fundamental_type
        {
            typedef data_t type;
        };
        template <has_value_type data_t> struct get_fundamental_type<data_t>
        {
            typedef typename data_t::value_type type;
        };
    }
    
    template <typename alias_t, typename mem_map_t, typename base_t, const array_centering ctr, typename grid_t>
    struct array_image_t
    {
        using alias_type       = alias_t;
        using fundamental_type = detail::get_fundamental_type<alias_t>::type;
        using index_type       = typename get_index_type<ctr>::array_type;
        using grid_type        = grid_t;
        using value_type       = fundamental_type;
        
        base_t      base;
        std::size_t csize;
        mem_map_t   map;
        
        constexpr static auto centering_type() { return ctr; }
        
        _sp_hybrid std::size_t size() const { return csize; }
        
        //Note: This has potential to become a bottleneck.
        template <typename... idxs_t>
        _sp_hybrid _finline_ const fundamental_type& operator() (const idxs_t&... idxs) const
        {
            return base[map.compute_offset(idxs...)];
        }
        
        template <typename... idxs_t>
        _sp_hybrid _finline_ fundamental_type& operator() (const idxs_t&... idxs)
        {
            return base[map.compute_offset(idxs...)];
        }
        
        _sp_hybrid _finline_ void set_elem(const index_type& idx, const alias_type& alias)
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int i = 0; i < alias.size(); ++i) (*this)(i, idx) = alias[i];
            }
            if constexpr (!ctrs::basic_array<alias_type>)
            {
                (*this)(idx) = alias;
            }
        }
        
        _sp_hybrid alias_type get_elem(const index_type& idx) const
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                alias_type output;
                for (int i = 0; i < output.size(); ++i) output[i] = (*this)(i, idx);
                return output;
            }
            else
            {
                return (*this)(idx);
            }
        }
        
        template <typename rhs_type>
        _sp_hybrid _finline_ void incr_elem(const index_type& idx, const rhs_type& rhs)
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int i = 0; i < alias_type::size(); ++i) (*this)(i, idx) += rhs[i];
            }
            else
            {
                (*this)(idx) += rhs;
            }
        }
        
        template <typename rhs_type>
        _sp_hybrid _finline_ void decr_elem(const index_type& idx, const rhs_type& rhs)
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int i = 0; i < alias_type::size(); ++i) (*this)(i, idx) -= rhs[i];
            }
            else
            {
                (*this)(idx) -= rhs;
            }
        }
        
        template <typename rhs_type>
        _sp_hybrid _finline_ void mult_elem(const index_type& idx, const rhs_type& rhs)
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int i = 0; i < alias_type::size(); ++i) (*this)(i, idx) *= rhs[i];
            }
            else
            {
                (*this)(idx) *= rhs;
            }
        }
        
        template <typename rhs_type>
        _sp_hybrid _finline_ void divi_elem(const index_type& idx, const rhs_type& rhs)
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int i = 0; i < alias_type::size(); ++i) (*this)(i, idx) /= rhs[i];
            }
            else
            {
                (*this)(idx) /= rhs;
            }
        }
    };
}