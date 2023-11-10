#pragma once

#include <string>
#include <sstream>

namespace spade::sym
{
    template <typename T, const T... ts>
    struct symbol_t
    {
        using designator_type = T;
        
        static std::string str()
        {
            std::stringstream oss;
            (oss << ... << ts);
            return oss.str();
        }
    };
}