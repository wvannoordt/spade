#pragma once

#include <iostream>
#include <concepts>
#include <type_traits>

namespace cvdf::dims
{
    template <class T> concept grid_array_dimension = requires(T t, std::size_t i)
    {
        t;
        T::rank();
        t.total_size();
        t.get_index_extent(i);
    };
    
    struct singleton_dim
    {
        static constexpr std::size_t rank(void) {return 0;}
        std::size_t total_size(void) const {return 1;}
        std::size_t get_index_extent(std::size_t i) {return 1;}
    };
    
    template <const std::size_t... s_dims_t> struct static_dims
    {
        static constexpr std::size_t rank(void) {return sizeof...(s_dims_t);};
        // std::size_t prod_r(void) const {return 1;}
        // template <const std::size_t p_dim_t> std::size_t prod_r(std::size_t d) const {return d;}
        template <typename... p_dims_t> std::size_t prod(p_dims_t... dims_in) const {return (... * dims_in);}
        std::size_t total_size(void) const {return prod(s_dims_t...);}
        template <typename p_dim_t> std::size_t get_elem_r(std::size_t i_cur, std::size_t i_target, p_dim_t dim_in)
        {
            if (i_cur == i_target) return dim_in;
            else return 0;
        }
        template <typename p_dim_t, typename... p_dims_t> std::size_t get_elem_r(std::size_t i_cur, std::size_t i_target, p_dim_t dim_in, p_dims_t... dims_in)
        {
            if (i_cur == i_target) return dim_in;
            else return get_elem_r(i_cur+1,i_target,dims_in...);
        }
        std::size_t get_index_extent(std::size_t i)
        {
            return get_elem_r(0,i,s_dims_t...);
        }
    };
    
    template <const std::size_t size_dim> struct dynamic_dims
    {
        static constexpr std::size_t rank(void) {return size_dim;}
        std::size_t dim_data[size_dim] = {0};
        dynamic_dims(){}
        template <typename sizes1_t> void set_r(std::size_t idx, sizes1_t size)
        {
            static_assert(std::is_integral<sizes1_t>::value, "dimensions must be integral");
            dim_data[idx]=size;
        }
        template <typename sizes1_t, typename... sizes_t> void set_r(std::size_t idx, sizes1_t size, sizes_t... sizes)
        {
            static_assert(std::is_integral<sizes1_t>::value, "dimensions must be integral");
            dim_data[idx]=size;
            set_r(idx+1, sizes...);
        }
        template <typename... sizes_t> dynamic_dims(sizes_t... sizes)
        {
            static_assert(sizeof...(sizes_t)==size_dim, "constructor requires number of arguments equal to rank");
            set_r(0, sizes...);
        }
        std::size_t get_index_extent(std::size_t i) {return dim_data[i];}
        std::size_t total_size(void) const
        {
            std::size_t output = 1;
            for (std::size_t i = 0; i < size_dim; i++) output*= dim_data[i];
            return output;
        }
    };
}