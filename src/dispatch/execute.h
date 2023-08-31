#pragma once

#include "core/except.h"
#include "core/ctrs.h"
#include "dispatch/device_type.h"
#include "dispatch/kernel_params.h"

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
    
    template <typename range_t, typename kernel_t>
    static void execute(const range_t& range, kernel_t& kernel)
    {
        using loop_index_type = typename range_t::index_type;
        const auto bound_box  = range.idx_bound_box();
        using bound_type = utils::remove_all<decltype(bound_box)>::type;
        if constexpr (device::is_cpu<typename range_t::device_t>)
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
            detail::k_gpu_dispatch_impl<<<ps.grid_size, ps.block_size>>>(ps, kernel);
            cudaDeviceSynchronize();
            auto er_code = cudaGetLastError();
            std::string er_str = std::string(cudaGetErrorString(er_code));
            if (er_code != cudaSuccess)
            {
                throw except::sp_exception("Error after kernel call: " + er_str);
            }
            //Do stuff
#endif
        }
    }
}