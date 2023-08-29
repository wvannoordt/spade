#pragma once

#include <vector>
#include <string>
#include "core/cuda_incl.h"
#include "dispatch/device_vector.h"

namespace spade::device
{
    template <typename data_t>
    struct shared_vector
    {
        using value_type = data_t;
        
        std::vector<data_t>   host_data;
        device_vector<data_t> devc_data;
        
        void push_back(const value_type& nval)
        {
            host_data.push_back(nval);
        }
        
        void resize(const std::size_t& n)
        {
            host_data.resize(n);
            devc_data.resize(n);
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