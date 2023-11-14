#pragma once

#include <tuple>
#include <string>
#include <sstream>

namespace spade::sym
{
    template <typename T> concept is_symbol = T::is_symbol_v;
    
    template <typename T, const T... ts>
    struct symbol_t
    {
        using designator_type = T;
        
        static constexpr bool is_symbol_v = true;
        
        static std::string str()
        {
            std::stringstream oss;
            (oss << ... << ts);
            return oss.str();
        }
    };
    
    template <typename T> concept is_vector = T::is_vector_v;
    
    template <is_symbol... vars_t>
    struct vector_t
    {
        using tuple_type = std::tuple<vars_t...>;
        static constexpr bool is_vector_v = true;
        vector_t (const vars_t&...){}
    };
}