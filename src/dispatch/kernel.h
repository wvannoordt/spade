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
        using exec_space_type = exec_space_t;
        using index_type      = typename outer_range_t::index_type;
        outer_range_t o_range;
        callable_t function;
        exec_space_t space;
        shmem_t shmem;
        
        _sp_hybrid std::size_t shmem_size() const
        {
            std::size_t output = shmem.bytes();
            if constexpr (!device::is_device<exec_space_t>) output += space.shmem_size();
            return output;
        }
        
        _sp_hybrid index_type compute_index(const ranges::d3_t& g_dim, const ranges::d3_t& b_idx) const
        {
            static_assert(outer_range_t::index_type::size() == 1, "only 1-dimensional inner range currently supported!");
            return b_idx.x;
        }
        
        _sp_hybrid bool valid(const index_type& idx) const
        {
            return o_range.bounds.contains(idx);
        }
    };
    
    template <typename range_t, typename func_t, typename dspace_t, typename shmem_t>
    inline auto make_kernel(const range_t& outer_range, func_t& func, const dspace_t& dspace, const shmem_t& shmem)
    {
        return kernel_t{outer_range, func, dspace, shmem};
    } 
}