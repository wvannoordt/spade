#pragma once

#include <tuple>
#include <string>
#include <sstream>

#include "core/utils.h"

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
        static constexpr bool is_vector_v = true;
        _sp_hybrid constexpr static int size() { return sizeof...(vars_t); }
        
        template <is_symbol sym_t>
        constexpr static int index_of = utils::index_of_pack<typename utils::remove_all<sym_t>::type, vars_t...>;
        
        template <is_symbol new_var_t>
        using extend = vector_t<vars_t..., new_var_t>;
        
        vector_t (const vars_t&...){}
    };
}