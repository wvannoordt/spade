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
        d3_t grid_dim {1,1,1};
        d3_t blck_dim {1,1,1};
        
        const auto i_div_up = [](int a, int b){ return (a % b != 0) ? (a / b + 1) : (a / b); };
        
        if constexpr (device::is_device<typename kernel_t::exec_space_type>)
        {
            //In this case, there is no "thread pool" object to deal with, so we just compute the
            //parameters as we normally would.
            if constexpr (kernel_t::index_type::size() == 1)
            {
                // 
                //  | + | + | + | + | + | + | + | + | + | + | + | + |
                
                grid_dim.x = i_div_up(kernel.o_range.bounds.volume(), 32);
                blck_dim.x = 32;
            }
            else
            {
                throw except::sp_exception("IMPLEMENTATION FOR RANGE NOT FOUND");
            }
        }
        else
        {
            const auto& thrds = kernel.space;
            //In this case, we are separating the kernel into an "outer" and "inner" range
            static_assert(kernel_t::exec_space_type::index_type::size() == 1, "only 1-dimensional inner range currently supported!");
            
            //This will need to be modified if we ever support multidimensional thread pools
            blck_dim.x = thrds.size();
            
            if constexpr (kernel_t::index_type::size() == 1)
            {
                grid_dim.x = kernel.o_range.bounds.volume();
            }
            else if constexpr (kernel_t::index_type::size() == 3)
            {
                grid_dim.x = kernel.o_range.bounds.size(0);
                grid_dim.y = kernel.o_range.bounds.size(1);
                grid_dim.z = kernel.o_range.bounds.size(2);
            }
            else
            {
                throw except::sp_exception("IMPLEMENTATION FOR RANGE NOT FOUND");
            }
        }
        
        return kernel_config_t{grid_dim, blck_dim, kernel.shmem_size()};
    }
}