#pragma once

#include <string>

#include "dispatch/device_type.h"
#include "dispatch/ranges/ranges.h"
#include "dispatch/ranges/kernel_config.h"
#include "dispatch/kernel_threads.h"
#include "dispatch/kernel.h"

namespace spade::dispatch::proto
{
    namespace detail
    {

#if (_sp_cuda)
        template <typename kernel_t>
        __global__ void k_gpu_impl(kernel_t kern)
        {
            //The size of this is a kernel launch parameter
            extern __shared__ char sh_mem_raw[];
            
            char* base        = (char*)sh_mem_raw;
            auto& thrds       = kern.space;
            thrds.set_shmem_base(base);
            char* shmem_start = base + thrds.shmem_size();
            auto& function    = kern.function;
            thrds.set_idx(threadIdx);
            function(0, thrds);
        }
#endif
    }
    
    template <typename range_t, typename kernel_t, typename comp_space_t, typename shmem_t = shmem::empty_t>
    inline void execute(const range_t& range, kernel_t& kernel_in, const comp_space_t& space, const shmem_t& shmem = shmem::empty_t())
    {
        if (range.is_empty()) return;
        //We compute the kernel launch parameters
        auto kernel = make_kernel(range, kernel_in, space, shmem);
        const auto config = ranges::compute_config(kernel);
        
        const auto gconf      = config.grid_dim.data;
        const auto bconf      = config.block_dim.data;
        const auto shmem_size = config.shmem_size;
        
#if (_sp_cuda)
        detail::k_gpu_impl<<<gconf, bconf, shmem_size, 0>>>(kernel);
        cudaDeviceSynchronize();
        auto er_code = cudaGetLastError();
        std::string er_str = std::string(cudaGetErrorString(er_code));
        if (er_code != cudaSuccess)
        {
            throw except::sp_exception("Error after kernel call: " + er_str);
        }
#endif
    }
    /*
    template <std::integral idx_t, const std::size_t ar_size, typename kernel_t, device::is_cpu device_t>
    inline void execute(const ranges::basic_range_t<idx_t, ar_size>& range, kernel_t& kernel, const device_t& device)
    {
        if (range.is_empty()) return;
    }
    */
}