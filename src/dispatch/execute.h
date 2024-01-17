#pragma once

#include "core/except.h"
#include "core/ctrs.h"
#include "dispatch/device_type.h"
#include "dispatch/kernel_params.h"
#include "dispatch/ranges/ranges.h"
#include "dispatch/ranges/kernel_config.h"
#include "dispatch/kernel_threads.h"
#include "dispatch/kernel.h"

namespace spade::dispatch
{
    namespace detail
    {
        template <typename index_t, typename bbx_type, const int idx, typename kernel_t>
        requires ((idx < 0) && ctrs::basic_array<index_t>)
        static void cpu_dispatch_impl(index_t& i, const bbx_type& bounds, kernel_t& kernel)
        {
            if constexpr (index_t::size() == 1) kernel(i[0]);
            else                                kernel(i);
        }
        
        template <typename index_t, typename bbx_type, const int idx, typename kernel_t>
        requires ((idx >= 0) && ctrs::basic_array<index_t>)
        static void cpu_dispatch_impl(index_t& i, const bbx_type& bounds, kernel_t& kernel)
        {
            for (i[idx] = bounds.min(idx); i[idx] < bounds.max(idx); ++i[idx])
            {
                cpu_dispatch_impl<index_t, bbx_type, idx-1, kernel_t>(i, bounds, kernel);
            }
        }
        
#if(_sp_cuda)
        
        template <typename launch_params_t, typename kern_t>
        __global__ void k_gpu_dispatch_impl(const launch_params_t mapping, kern_t kern)
        {
            const auto i = mapping.get_index(threadIdx, blockIdx, blockDim, gridDim);
            if (mapping.is_valid(i))
            {
                if constexpr (launch_params_t::range_type::index_type::size() == 1)
                {
                    kern(i[0]);
                }
                else
                {
                    kern(i);
                }
            }
        }

#endif
    }
    
    template <typename range_t, typename kernel_t, typename device_t>
    static void execute(const range_t& range, kernel_t& kernel, const device_t& device)
    {
        using loop_index_type = typename range_t::index_type;
        const auto bound_box  = [&]()
        {
            if constexpr (requires {range.idx_bound_box();}) return range.idx_bound_box();
            else return range;
        }();
        
        if (bound_box.volume() == 0) return;
        
        using bound_type = utils::remove_all<decltype(bound_box)>::type;
        if constexpr (device::is_cpu<device_t>)
        {
            //CPU dispatch
            loop_index_type i;
            detail::cpu_dispatch_impl<
                loop_index_type,
                bound_type,
                loop_index_type::size()-1,
                kernel_t>(i, bound_box, kernel);
        }
        else
        {
            //GPU dispatch
            const auto ps = get_launch_params(range);
#if(_sp_cuda)
            ///Notes on kernel dispatch
            // k_kernel<<<grid_size, block_size, shmem_size, istream>>>(...);
            // grid_size  = size of CUDA grid in blocks   (grid_size.x,  grid_size.y,  grid_size.z)
            // block_size = size of CUDA block in threads (block_size.x, block_size.y, block_size.z)
            // shmem_size = size of shared memory allocated (in bytes) allocated per block
            // istream    = ID of the cuda stream that is dispatched
            
            //generalize this later
            // print(ps.grid_size.x, ps.grid_size.y, ps.grid_size.z);
            // print(ps.block_size.x, ps.block_size.y, ps.block_size.z);
            // std::cin.get();
            detail::k_gpu_dispatch_impl<<<ps.grid_size, ps.block_size>>>(ps, kernel);
            cudaDeviceSynchronize();
            auto er_code = cudaGetLastError();
            std::string er_str = std::string(cudaGetErrorString(er_code));
            if (er_code != cudaSuccess)
            {
                throw except::sp_exception("Error after kernel call: " + er_str);
            }
#endif
        }
    }
    
    template <typename range_t, typename kernel_t>
    requires (device::is_cpu<typename range_t::device_t> || device::is_gpu<typename range_t::device_t>)
    static void execute(const range_t& range, kernel_t& kernel)
    {
        execute(range, kernel, range.device());
    }
    
    
    
    //generalized ranges and shared memory
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
            const auto index  = kern.compute_index(gridDim, blockIdx);
            
            constexpr bool has_shmem = !(kernel_t::shared_type::is_empty());
            if constexpr (has_shmem)
            {
                kern.shmem.bind_ptr(shmem_start);
            }
            if (thrds.valid() && kern.valid(index))
            {
                if constexpr (has_shmem)
                {
                    if constexpr (kernel_t::index_type::size() == 1) function(index[0], thrds, kern.shmem);
                    else function(index, thrds, kern.shmem);
                }
                else
                {
                    if constexpr (kernel_t::index_type::size() == 1) function(index[0], thrds);
                    else function(index, thrds);
                }
            }
        }
#endif
    }
    
    template <std::integral idx_t, const std::size_t ar_size, typename kernel_t, typename comp_space_t, typename shmem_t = shmem::empty_t>
    inline void execute(const ranges::basic_range_t<idx_t, ar_size>& range, kernel_t& kernel_in, const comp_space_t& space, const shmem_t& shmem = shmem::empty_t())
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
}