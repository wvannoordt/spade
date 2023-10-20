#pragma once

#include <vector>
#include <string>
#include "core/cuda_incl.h"
#include "dispatch/device_vector.h"

namespace spade::device
{
    // Note that this concept should succeed for std::vector<...>
    // and shared_vector<...>, BUT NOT device_vector<...> (you can't push_back on device)
    template <typename T> concept sizeable_vector = requires (T t)
    {
        t.resize(1UL);
        t.push_back(typename T::value_type());
        t.size();
        t[0UL];
    };
    
    template <typename data_t>
    struct shared_vector
    {
        using value_type = data_t;
        
        std::vector<data_t>   host_data;
        device_vector<data_t> devc_data; // we will make this an array later when we do multiple gpus
        
        template <device::is_device device_t>
        auto& data(const device_t&)
        {
            if constexpr (device::is_cpu<device_t>) return host_data;
            else                                    return devc_data;
        }
        
        template <device::is_device device_t>
        const auto& data(const device_t&) const
        {
            if constexpr (device::is_cpu<device_t>) return host_data;
            else                                    return devc_data;
        }
        
        void push_back(const value_type& nval)
        {
            host_data.push_back(nval);
        }
        
        void resize(const std::size_t& n)
        {
            host_data.resize(n);
            devc_data.resize(n);
        }
        
        auto size() const
        {
            return host_data.size();
        }
        
        auto begin() const
        {
            return host_data.begin();
        }
        
        auto end() const
        {
            return host_data.end();
        }
        
        void clear()
        {
            this->resize(0);
        }
        
        void transfer()
        {
#if(_sp_cuda)
            devc_data.resize(host_data.size());
            std::string err_string = "no error";
            void*       dest = (void*)(devc_data.raw);
            const void* src  = (const void*)(&host_data[0]);
            auto er_code = cudaMemcpy(dest, src, host_data.size() * sizeof(data_t), cudaMemcpyHostToDevice);
            if (er_code == cudaSuccess) return;
            err_string = std::string(cudaGetErrorString(er_code));
            throw mem_exception("unable to transfer to device: " + err_string);
#endif
        }
        
        void itransfer()
        {
#if(_sp_cuda)
            host_data.resize(devc_data.size());
            std::string err_string = "no error";
            void*       dest = (void*)(host_data.data());
            const void* src  = (const void*)(&devc_data[0]);
            auto er_code = cudaMemcpy(dest, src, host_data.size() * sizeof(data_t), cudaMemcpyDeviceToHost);
            if (er_code == cudaSuccess) return;
            err_string = std::string(cudaGetErrorString(er_code));
            throw mem_exception("unable to transfer to device: " + err_string);
#endif
        }
        
        data_t&       operator[] (const std::size_t& idx)       { return host_data[idx]; }
        const data_t& operator[] (const std::size_t& idx) const { return host_data[idx]; }
    };
}