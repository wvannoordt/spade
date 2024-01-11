#pragma once

#include "dispatch/shmem.h"

namespace spade::dispatch
{
    // Note: this needs to be value-copyable
    //
    // callable_t:   the actual kernel code, consisting of the lambda
    // exec_space_t: either cpu/gpu OR kernel_threads_t<gpu/cpu>.
    //
    template <typename outer_range_t, typename callable_t, typename exec_space_t, typename shmem_t>
    struct kernel_t
    {
        outer_range_t o_range;
        callable_t function;
        exec_space_t space;
        shmem_t shmem;
        
        std::size_t shmem_size() const
        {
            std::size_t output = shmem.bytes();
            if constexpr (!device::is_device<exec_space_t>) output += space.shmem_size();
            return output;
        }
    };
    
    template <typename range_t, typename func_t, typename dspace_t, typename shmem_t>
    inline auto make_kernel(const range_t& outer_range, func_t& func, const dspace_t& dspace, const shmem_t& shmem)
    {
        return kernel_t{outer_range, func, dspace, shmem};
    } 
}