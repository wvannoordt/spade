#pragma once

#include <concepts>

#include "core/ctrs.h"
//TODO

namespace spade::block_config
{
    template <typename T> concept block_config = requires(T t)
    {
        t;
    };
    
    template <typename coord_val_t, const std::size_t dim> struct cartesian_blocks_t
    {
        cartesian_blocks_t(const ctrs::array<int, dim>& num_blocks_in, const spade::bound_box_t<coord_val_t, dim>& bounds_in)
        {
            
        }
    };
}