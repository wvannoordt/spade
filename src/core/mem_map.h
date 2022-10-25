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
    
    template <const default_index_t i0, const default_index_t i1>
    static std::ostream & operator<<(std::ostream & os, const static_dim_t<i0, i1>& dim)
    {
        os << "static dim: [" << dim.first() << " --> " << dim.last() << "]";
        return os;
    }
    
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
        _finline_ offset_type size() const {return last()-first();}
        _finline_ offset_type offset(const index_type& i) const {return i-first();}
        
    };
    
    template <typename idx_t>
    static std::ostream & operator<<(std::ostream & os, const dynamic_dim_t<idx_t>& dim)
    {
        os << "dynamic dim: [" << dim.first() << " --> " << dim.last() << "]";
        return os;
    }
    
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
    
    static std::ostream & operator<<(std::ostream & os, const singleton_map_t& map)
    {
        os << "[singleton map]\n";
        os << "[/singleton map]";
        return os;
    }
    
    template <is_dimension... dims_t> struct regular_map_t
    {
        using index_type  = default_index_t; //TODO: change this to be promotion from the index types of the dimensions (later)
        using offset_type = default_offset_t;
        
        constexpr static index_type rank() {return sizeof...(dims_t);}
        
        using identifier_type = ctrs::array<index_type, rank()>;
        
        ctrs::array<offset_type, rank()> coeffs;
        aliases::tuple<dims_t...> dims;
        
        std::size_t map_size;
        
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
                coeffs[i.value] = coeffs[i.value-1]*(std::get<i.value-1>(dims).size());
            });
            map_size = 1;
            static_for<0,rank()>([&](const auto& i) -> void
            {
                map_size *= std::get<i.value>(dims).size();
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
        
        offset_type size() const
        {
            return map_size;
        }
    };
    
    template <typename... dims_t>
    static std::ostream & operator<<(std::ostream & os, const regular_map_t<dims_t...>& map)
    {
        os << "[regular map]\n";
        static_for<0,sizeof...(dims_t)>([&](const auto& i) -> void
        {
            os << std::get<i.value>(map.dims);
            os << "\n";
        });
        os << "[/regular map]";
        return os;
    }
    
    
    template <is_map... maps_t> struct composite_map_t
    {
        using index_type  = default_index_t;
        using offset_type = default_offset_t;
        
        static constexpr index_type rank() {return sizeof...(maps_t);}
        aliases::tuple<maps_t...> maps;
        ctrs::array<offset_type, rank()> coeffs;
        std::size_t map_size;
        
        composite_map_t()
        {
            compute_coeffs();
        }
        
        composite_map_t(maps_t... maps_in) //indentionally not const ref!
        {
            maps = std::make_tuple(maps_in...);
            compute_coeffs();
        }
        
        void compute_coeffs()
        {
            coeffs[0] = 1;
            static_for<1,rank()>([&](const auto& i) -> void
            {
                coeffs[i.value] = std::get<i.value-1>(maps).size()*coeffs[i.value-1];
            });
            map_size = 1;
            static_for<0, rank()>([&](const auto& i) -> void
            {
                map_size *= std::get<i.value>(maps).size();
            });
        }
        
        offset_type size() const
        {
            return map_size;
        }
        
        //increment using an array index with no further arguments
        template <const int cur_map, const int counter, ctrs::basic_array arg_t>
        requires (
            (arg_t::size() > 0) &&
            (utils::get_pack_type<cur_map, maps_t...>::type::rank() == arg_t::size()) &&
            (utils::get_pack_type<cur_map, maps_t...>::type::rank()>0))
        _finline_ void recurse_incr_offset(
            offset_type& cur_offset,
            const arg_t& arg) const
        {
            cur_offset += coeffs[cur_map]*std::get<cur_map>(maps).offset(arg);
        }
        
        //increment using an array index
        template <const int cur_map, const int counter, typename arg_t, typename... args_t>
        requires (
            (arg_t::size() > 0)
            && (utils::get_pack_type<cur_map, maps_t...>::type::rank() == arg_t::size()) &&
            (utils::get_pack_type<cur_map, maps_t...>::type::rank()>0))
        _finline_ void recurse_incr_offset(
            offset_type& cur_offset,
            const arg_t& arg,
            const args_t&... args) const
        {
            cur_offset += coeffs[cur_map]*std::get<cur_map>(maps).offset(arg);
            recurse_incr_offset<cur_map+1, 0, args_t...>(cur_offset, args...);
        }
        
        //Pull the first argument into a newly-created array
        template <const int cur_map, const int counter, typename arg_t, typename... args_t>
        requires (
            !ctrs::basic_array<arg_t>
            && counter==0
            && utils::get_pack_type<cur_map, maps_t...>::type::rank()>0 &&
            (utils::get_pack_type<cur_map, maps_t...>::type::rank()>0))
        _finline_ void recurse_incr_offset(
            offset_type& cur_offset,
            const arg_t& arg,
            const args_t&... args) const
        {
            using idx_t = utils::get_pack_type<cur_map, maps_t...>::type::identifier_type;
            idx_t ar;
            ar[0] = arg;
            recurse_incr_offset<cur_map, counter+1, args_t...>(cur_offset, ar, args...);
        }
        
        //Pull an argument into an accumulating array
        template <const int cur_map, const int counter, typename arg_t, typename... args_t>
        requires (
            !ctrs::basic_array<arg_t> &&
            counter>0 &&
            ((utils::get_pack_type<cur_map, maps_t...>::type::rank())>counter) &&
            (utils::get_pack_type<cur_map, maps_t...>::type::rank()>0))
        _finline_ void recurse_incr_offset(
            offset_type& cur_offset,
            typename utils::get_pack_type<cur_map, maps_t...>::type::identifier_type& ar,
            const arg_t& arg,
            const args_t&... args) const
        {
            ar[counter] = arg;
            recurse_incr_offset<cur_map, counter+1, args_t...>(cur_offset, ar, args...);
        }
        
        //Send a full array to the indexing overload
        // template <const int cur_map, const int counter, typename arg_t, typename... args_t>
        template <const int cur_map, const int counter, typename... args_t>
        requires (
            (counter >= utils::get_pack_type<cur_map, maps_t...>::type::rank()) &&
            (utils::get_pack_type<cur_map, maps_t...>::type::rank()>0))
        _finline_ void recurse_incr_offset(
            offset_type& cur_offset,
            typename utils::get_pack_type<cur_map, maps_t...>::type::identifier_type& ar,
            // const arg_t& arg,
            const args_t&... args) const
        {
            recurse_incr_offset<
                cur_map,
                0,
                typename utils::get_pack_type<cur_map, maps_t...>::type::identifier_type,
                args_t...>(cur_offset, ar, args...);
        }
        
        //increment past singleton maps
        template <const int cur_map, const int counter, typename... args_t>
        requires (utils::get_pack_type<cur_map, maps_t...>::type::rank() == 0)
        _finline_ void recurse_incr_offset(
            offset_type& cur_offset,
            const args_t&... args) const
        {
            recurse_incr_offset<
                cur_map+1,
                counter,
                args_t...>(cur_offset, args...);
        }
        
        template <typename... idxs_t> 
        _finline_ offset_type offset(const idxs_t&... idxs) const
        {
            offset_type output = 0;
            recurse_incr_offset<0,0,idxs_t...>(output, idxs...);
            return output;
        }
    };
    
    template <typename... maps_t>
    static std::ostream & operator<<(std::ostream & os, const composite_map_t<maps_t...>& map)
    {
        os << "[composite map]\n";
        static_for<0,sizeof...(maps_t)>([&](const auto& i) -> void
        {
            os << std::get<i.value>(map.maps);
            os << "\n";
        });
        os << "[/composite map]";
        return os;
    }
    
    template <typename... inputs_t> struct something_to_check_validity_of_index_inputs_here
    {
        static constexpr bool value = true;
    };
}