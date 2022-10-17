#pragma once

#include <type_traits>
#include <concepts>

#include "core/aliases.h" 

namespace spade::grid
{
    template <typename T> concept is_dimension = requires(T t, const typename T::index_type& i)
    {
        //The starting index dimension (inclusive)
        t.first();
        
        //The ending index dimension (exclusive)
        t.last();
    };
    
    template <typename T> concept is_map = requires(T t)
    {
        //The total size of the iteration space
        t.size();
    };
    
    using default_index_t = int;
    using default_offset_t = std::size_t;
    
    template <typename... maps_t> struct composite_map_t
    {
        using index_type = int;
        using offset_type = std::size_t;
        
        composite_map_t(){}
        
        template <typename... idxs_t> offset_type compute_offset(const idxs_t&... idxs)
        {
            return 0;
        }
    };
    
    template <const default_index_t i0, const default_index_t i1> struct static_dim_t
    {
        using index_type = default_index_t;
        using offset_type = default_offset_t;
        
        static_dim_t(){}
        
        constexpr index_type first() const {return i0;}
        constexpr index_type last()  const {return i1;}
    };
    
    template <typename idx_t>
    requires (std::is_integral_v<idx_t>)
    struct dynamic_dim_t
    {
        using index_type = idx_t;
        using offset_type = default_offset_t;
        
        idx_t i0, i1;
        
        dynamic_dim_t(){}
        
        index_type first() const {return i0;}
        index_type last()  const {return i1;}
    };
    
    template <is_dimension dim_t> static constexpr auto dim_size(const dim_t& dim)
    {
        return typename dim_t::offset_type(dim.first()-dim.last());
    }
    
    template <is_dimension... dims_t> struct regular_map_t
    {
        aliases::tuple<dims_t...> dims;
        regular_map_t(){}
        regular_map_t(dims_t... dims_in) //note that these are deliberately by value and not const ref!!!!!
        {
            dims = std::make_tuple(dims_in...);
        }
    };
    
    template <typename... inputs_t> struct something_to_check_validity_of_index_inputs_here
    {
        static constexpr bool value = true;
    };
}