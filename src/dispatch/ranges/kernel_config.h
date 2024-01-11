#pragma once

#include "core/cuda_incl.h"

namespace spade::dispatch::ranges
{
    struct dim3_wrapper_t
    {
        int x, y, z;
    };
    
#if (_sp_cuda)
    using d3_t = dim3;
#else
    using d3_t = dim3_wrapper_t;
#endif
    
    struct outer_config_t { d3_t data; };
    struct inner_config_t { d3_t data; };
    
    // <<<grid, block, shmem_sz, strm>>>
    //Later, we need to use a cuda stream too
    struct kernel_config_t
    {
        outer_config_t grid_dim;
        inner_config_t block_dim;
        std::size_t shmem_size;
    };
    
    template <typename kernel_t>
    auto compute_config(const kernel_t& kernel)
    {
        d3_t out0 {1,1,1};
        d3_t out1 {1,1,1};
        return kernel_config_t{out0, out1, kernel.shmem_size()};
    }
}