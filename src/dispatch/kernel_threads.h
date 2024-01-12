#pragma once

#include "core/vec_image.h"

#include "dispatch/device_type.h"
#include "dispatch/ranges/kernel_config.h"

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
            mutable volatile char* shmem_base = nullptr;
            
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
            _sp_hybrid std::size_t shmem_size() const { return size()*sizeof(double); }
            
            _sp_hybrid void set_shmem_base(volatile char* bin) { shmem_base = bin; }
            
            _sp_hybrid device_t device() const { return k_device; }
            
            _sp_hybrid void set_idx(const ranges::d3_t& b_idx)
            {
                static_assert(index_type::size() == 1, "only 1-dimensional inner range currently supported!");
                k_index = b_idx.x;
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
                    static_assert(index_type::size() == 1, "only 1-dimensional inner range currently supported!");
                    for (index_type i{irange.bounds.min(0)}; i<irange.bounds.max(0); ++i)
                    {
                        func(this->reduce_idx(i));
                    }
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
            
            template <const k_thread_reduce_op oper, typename data_t, typename func_t>
            _sp_hybrid inline auto threads_reduce(const reduction_t<oper, data_t>&, const func_t& func) const
            {
                using output_t = data_t;
                if constexpr (is_gpu_impl)
                {
                    utils::vec_image_t buf {(output_t*)shmem_base, this->size()};
                    static_assert(index_type::size() == 1, "only 1-dimensional inner range currently supported!");
                    buf[this->reduce_idx(k_index)] = func(this->reduce_idx(k_index));
                    
                    //This is probably a bit slow and janky! Note: we should do this with interleaving
                    //if it becomes a bottleneck
                    static_assert(index_type::size() == 1, "only 1-dimensional inner range currently supported!");
                    const auto root = irange.bounds.min(0);
                    if (k_index == root)
                    {
                        for (typename index_type::value_type i = root; i < irange.bounds.max(0); ++i)
                        {
                            if constexpr (oper == k_reduce_sum) buf[root] += (i!=root)*buf[i];
                            if constexpr (oper == k_reduce_max) buf[root] = utils::min(buf[root], buf[i]);
                            if constexpr (oper == k_reduce_min) buf[root] = utils::max(buf[root], buf[i]);
                        }
                    }
                    this->sync();
                    return buf[root];
                }
                else
                {
                    
                }
                
                return output_t();
            }
            
            template <typename func_t>
            _sp_hybrid inline auto sum(const func_t& func) const
            {
                using data_t = decltype(func(this->reduce_idx(index_type())));
                return this->threads_reduce(reduction_t<k_reduce_sum, data_t>(), func);
            }
            
            template <typename func_t>
            _sp_hybrid inline auto min(const func_t& func) const
            {
                using data_t = decltype(func(this->reduce_idx(index_type())));
                return this->threads_reduce(reduction_t<k_reduce_min, data_t>(), func);
            }
            
            template <typename func_t>
            _sp_hybrid inline auto max(const func_t& func) const
            {
                using data_t = decltype(func(this->reduce_idx(index_type())));
                return this->threads_reduce(reduction_t<k_reduce_max, data_t>(), func);
            }
    };
}