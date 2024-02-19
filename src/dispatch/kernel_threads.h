#pragma once

#include "core/vec_image.h"

#include "dispatch/device_type.h"
#include "dispatch/ranges/kernel_config.h"
#include "dispatch/range_loop.h"

namespace spade::dispatch
{
    enum k_thread_reduce_op
    {
        k_reduce_min,
        k_reduce_max,
        k_reduce_sum
    };
    
    template <const k_thread_reduce_op, typename data_t> struct reduction_t {};
    
    template <typename irange_t, device::is_device device_t>
    class kernel_threads_t
    {
        public:
            using index_type = irange_t::index_type;
        
        constexpr static bool is_gpu_impl = device::is_gpu<device_t>;
        
        private:
            //This might be absolutely terrible
            index_type k_index;
            
            _sp_device void i_sncthr() const
            {
#if (_sp_cuda)
                __syncthreads();
#endif    
            }
            
            _sp_hybrid inline auto reduce_idx(const index_type& idx) const
            {
                if constexpr (index_type::size() == 1)
                {
                    return idx[0];
                }
                else
                {
                    return idx;
                }
            }
        
        public:
            irange_t irange;
            device_t k_device;
            kernel_threads_t(const irange_t& irange_in, const device_t& device_in) : irange{irange_in}, k_device{device_in} {}
            
            _sp_hybrid std::size_t size() const { return irange.bounds.volume(); }
            
            // This is a goofy idea so far...
            // _sp_hybrid std::size_t shmem_size() const { return size()*sizeof(double); }
            
            // _sp_hybrid void set_shmem_base(volatile char* bin) const { shmem_base = bin; }
            
            _sp_hybrid device_t device() const { return k_device; }
            
            _sp_hybrid bool isroot() const
            {
                if constexpr (is_gpu_impl)
                {
                    k_index == 0;
                }
                else
                {
                    return true;
                }
            }
            
            _sp_hybrid void set_idx(const ranges::d3_t& b_idx)
            {
                static_assert((index_type::size() == 1) || (index_type::size() == 2) || (index_type::size() == 3), "only 1/2/3-dimensional inner range currently supported!");
                if constexpr (index_type::size() == 1)
                {
                    k_index = b_idx.x;
                }
                if constexpr (index_type::size() == 2)
                {
                    k_index[0] = b_idx.x;
                    k_index[1] = b_idx.y;
                }
                else if constexpr (index_type::size() == 3)
                {
                    k_index[0] = b_idx.x;
                    k_index[1] = b_idx.y;
                    k_index[2] = b_idx.z;
                }
            }
            
            template <typename func_t>
            _sp_hybrid void exec(const func_t& func) const
            {
                if constexpr (is_gpu_impl)
                {
                    func(this->reduce_idx(k_index));
                }
                else
                {
                    //TODO:: replace later with the cpu for-loop for multidimensional range
                    dispatch::range_loop(irange, func);
                    // for (typename index_type::value_type i{irange.bounds.min(0)}; i<irange.bounds.max(0); ++i)
                    // {
                    //     func(this->reduce_idx(i));
                    // }
                }
            }
            
            _sp_hybrid void sync() const
            {
                if constexpr (is_gpu_impl)
                {
                    i_sncthr();
                }
            }
            
            _sp_hybrid bool valid() const
            {
                return irange.bounds.contains(k_index);
            }
    };
}