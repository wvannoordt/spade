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
        
        if constexpr (device::is_device<typename kernel_t::exec_space_type>)
        {
            //In this case, there is no "thread pool" object to deal with, so we just compute the
            //parameters as we normally would.
            if constexpr (kernel_t::index_type::size() == 1)
            {
                // 
                //  | + | + | + | + | + | + | + | + | + | + | + | + |
                using intgr_t = decltype(kernel.o_range.bounds.volume());
                grid_dim.x = utils::i_div_up(kernel.o_range.bounds.volume(), intgr_t(32));
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
            static_assert(
                (kernel_t::exec_space_type::index_type::size() == 1) ||
                (kernel_t::exec_space_type::index_type::size() == 2) ||
                (kernel_t::exec_space_type::index_type::size() == 3),
                "only 1 or 2 or 3-dimensional inner range currently supported!");
            
            //This will need to be modified if we ever support multidimensional thread pools
            if constexpr (kernel_t::exec_space_type::index_type::size() == 1)
            {
                blck_dim.x = thrds.size();
            }
            else if constexpr (kernel_t::exec_space_type::index_type::size() == 2)
            {
                blck_dim.x = thrds.irange.bounds.size(0);
                blck_dim.y = thrds.irange.bounds.size(1);
            }
            else if constexpr (kernel_t::exec_space_type::index_type::size() == 3)
            {
                blck_dim.x = thrds.irange.bounds.size(0);
                blck_dim.y = thrds.irange.bounds.size(1);
                blck_dim.z = thrds.irange.bounds.size(2);
            }
            
            if constexpr (kernel_t::index_type::size() == 1)
            {
                grid_dim.x = kernel.o_range.bounds.volume();
            }
            else if constexpr (kernel_t::index_type::size() == 2)
            {
                grid_dim.x = kernel.o_range.bounds.size(0);
                grid_dim.y = kernel.o_range.bounds.size(1);
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