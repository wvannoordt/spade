#pragma once

#include <concepts>

#include "core/ctrs.h"
//TODO

namespace spade::block_config
{
    template <typename T> concept block_configuration = requires(T t)
    {
        typename T::node_type;
    };
    
    template <typename coord_val_t, const std::size_t grid_dim> struct cartesian_blocks_t
    {
        ctrs::array<int, grid_dim> num_blocks;
        bound_box_t<coord_val_t, grid_dim> bounds;
        
        cartesian_blocks_t(const ctrs::array<int, grid_dim>& num_blocks_in, const bound_box_t<coord_val_t, grid_dim>& bounds_in)
        {
            num_blocks = num_blocks_in;
            bounds     = bounds_in;
        }
        
    };
}