#pragma once

namespace spade::time_integration
{
    struct rk2_t
    {
        static constexpr std::size_t var_size() {return 1;}
        static constexpr std::size_t rhs_size() {return 2;}
    };
    
    struct rk4_t
    {
        static constexpr std::size_t var_size() {return 3;}
        static constexpr std::size_t rhs_size() {return 4;}
    };
}