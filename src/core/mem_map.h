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
        
        t.size();
        t.offset(i);
    };
    
    template <typename T> concept is_map = requires(T t)
    {
        //The total size of the iteration space
        t.size();
        
        //number of indices required to produce an offset
        t.rank();
    };
    
    using default_index_t = int;
    using default_offset_t = std::size_t;
    
    template <const default_index_t i0, const default_index_t i1> struct static_dim_t
    {
        using index_type = default_index_t;
        using offset_type = default_offset_t;
        
        static_dim_t(){}
        
        constexpr index_type first() const {return i0;}
        constexpr index_type last()  const {return i1;}
        constexpr offset_type size() const {return last()-first();}
        constexpr offset_type offset(const index_type& i) const {return i-first();}
    };
    
    template <typename idx_t>
    requires (std::is_integral_v<idx_t>)
    struct dynamic_dim_t
    {
        using index_type = idx_t;
        using offset_type = default_offset_t;
        
        idx_t i0, i1;
        
        dynamic_dim_t(){}
        dynamic_dim_t(const idx_t& i0_in, const idx_t& i1_in){ i0 = i0_in; i1 = i1_in; }
        index_type first() const {return i0;}
        index_type last()  const {return i1;}
        offset_type size() const {return last()-first();}
        offset_type offset(const index_type& i) const {return i-first();}
        
    };
    
    struct singleton_map_t
    {
        using index_type  = default_index_t;
        using offset_type = default_offset_t;
        
        static constexpr offset_type size() {return 1;}
        static constexpr index_type  rank() {return 0;}
        
        using identifier_type = ctrs::array<index_type, 0>;
        
        _finline_ offset_type offset(const identifier_type& idx) const
        {
            return 0;
        }
    };
    
    template <is_dimension... dims_t> struct regular_map_t
    {
        using index_type  = default_index_t; //TODO: change this to be promotion from the index types of the dimensions (later)
        using offset_type = default_offset_t;
        
        constexpr static index_type rank() {return sizeof...(dims_t);}
        
        using identifier_type = ctrs::array<index_type, rank()>;
        
        ctrs::array<offset_type, rank()> coeffs;
        aliases::tuple<dims_t...> dims;
        
        regular_map_t()
        {
            //maybe do this only when the sizes are constant-eval
            compute_coeffs();
        }
        
        regular_map_t(dims_t... dims_in) //note that these are deliberately by value and not const ref!!!!!
        {
            dims = std::make_tuple(dims_in...);
            compute_coeffs();
        }
        
        void compute_coeffs()
        {
            coeffs[0] = 1;
            static_for<1,rank()>([&](const auto& i) -> void
            {
                coeffs[i.value] = coeffs[i.value-1]*(std::get<i.value>(dims).size());
            });
        }
        
        _finline_ offset_type offset(const identifier_type& idx) const
        {
            offset_type output = offset_type(0);
            static_for<0,rank()>([&](const auto& i) -> void
            {
                output += std::get<i.value>(dims).offset(idx[i.value])*coeffs[i.value];
            });
            return output;
        }
    };
    
    
    template <typename... maps_t> struct composite_map_t
    {
        using index_type  = default_index_t;
        using offset_type = default_offset_t;
        
        aliases::tuple<maps_t...> maps;
        composite_map_t(){}
        composite_map_t(maps_t... maps_in) //indentionally not const ref!
        {
            maps = std::make_tuple(maps_in...);
        }
        
        template <const int cur_map, const int s_count, typename arg_t, typename... args_t>
        requires ((s_count == 0) && !ctrs::basic_array<arg_t>)
        void recurse_incr_offset(offset_type& cur_offset, const arg_t& arg, const args_t&... args) const
        {
            using map_t = utils::get_pack_type<cur_map, maps_t...>::type;
            using ar_t  = map_t::identifier_type;
            ar_t ar;
            //START HERE
            return 0;
        }
        
        template <typename... idxs_t> offset_type offset(const idxs_t&... idxs) const
        {
            offset_type output = 0;
            return recurse_incr_offset<0,0,idxs_t...>(output, idxs...);
        }
    };
    
    template <typename... inputs_t> struct something_to_check_validity_of_index_inputs_here
    {
        static constexpr bool value = true;
    };
}