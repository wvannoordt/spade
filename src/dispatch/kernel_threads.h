#pragma once

#include "dispatch/device_type.h"
#include "dispatch/ranges/kernel_config.h"

namespace spade::dispatch
{
    template <typename irange_t, device::is_device device_t>
    class kernel_threads_t
    {
        using index_type = irange_t::index_type;
        
        constexpr static bool is_gpu_impl = device::is_gpu<device_t>;
        
        private:
            //This might be absolutely terrible
            index_type k_index;
            volatile char* shmem_base = nullptr;
        
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
                    func(k_index);
                }
                else
                {
                    for (index_type i{irange.bounds.min(0)}; i<irange.bounds.max(0); ++i)
                    {
                        func(i);
                    }
                }
            }
    };
}