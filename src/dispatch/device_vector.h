#pragma once

#include <vector>
#include "dispatch/device_type.h"
#include "dispatch/device_allocator.h"

namespace spade::device
{
    //todo: improve this implementation
    template <typename data_t, typename allocator_t = device_allocator_t<data_t>>
    struct device_vector
    {
        using value_type = data_t;
        allocator_t allocator;
        data_t* raw = nullptr;
        std::size_t c_size = 0;
        
        device_vector(){}
        device_vector(const std::size_t n) 
        {
            this->resize(n);
        }
        
        device_vector(const device_vector& rhs) 
        {
            if (rhs.size() != 0)
            {
                this->resize(rhs.size());
#if (_sp_cuda)
                cudaMemcpy(this->raw, rhs.raw, c_size*sizeof(value_type), cudaMemcpyDeviceToDevice);
#endif
            }
        }
        
        operator std::vector<data_t>() const
        {
            std::vector<data_t> out(c_size);
            data_t* base = out.data();
#if (_sp_cuda)
            cudaMemcpy(base, raw, c_size*sizeof(value_type), cudaMemcpyDeviceToHost);
#endif
            return out;
        }
        
        device_vector& operator = (const device_vector& rhs)
        {
            if (rhs.size() == 0) return *this;
            this->resize(rhs.size());
#if (_sp_cuda)
            if (rhs.size() > 0) cudaMemcpy(this->raw, rhs.raw, c_size*sizeof(value_type), cudaMemcpyDeviceToDevice);
#endif
            return *this;
        }
    
        device_vector& operator = (const std::vector<data_t>& rhs)
        {
            if (rhs.size() == 0) return *this;
            this->resize(rhs.size());
#if (_sp_cuda)
            if (rhs.size() > 0) cudaMemcpy(this->raw, rhs.data(), c_size*sizeof(value_type), cudaMemcpyHostToDevice);
#endif
            return *this;
        }
        
        device_vector(device_vector&& rhs)
        {
            clear();
            if (rhs.size() != 0)
            {
                this->raw    = rhs.raw;
                this->c_size = rhs.c_size;
                rhs.raw      = nullptr;
                rhs.c_size   = 0;
            }
        }
        
        device_vector& operator = (device_vector&& rhs)
        {
            clear();
            if (rhs.size() != 0)
            {
                this->raw    = rhs.raw;
                this->c_size = rhs.c_size;
                rhs.raw      = nullptr;
                rhs.c_size   = 0;
            }
            return *this;
        }
        
        std::size_t size() const {return c_size;}
        
        void resize(const std::size_t& n)
        {
            if (n == c_size) return;
#if (_sp_cuda)
            if (raw != nullptr) allocator.deallocate(raw, c_size);
            c_size = n;
            raw = allocator.allocate(n);
#endif
        }
        
        void clear()
        {
            this->resize(0);
        }
        
        ~device_vector()
        {
            if (raw != nullptr) allocator.deallocate(raw, c_size);
        }
        
        _sp_hybrid data_t& operator[]       (const std::size_t& idx)       { return raw[idx]; }
        _sp_hybrid const data_t& operator[] (const std::size_t& idx) const { return raw[idx]; }
    };
    
    template <typename data_t, typename device_t>
    using auto_vector = std::conditional<is_cpu<device_t>, std::vector<data_t>, device_vector<data_t>>::type;
    
    template <typename data_t, typename allocator_t>
    inline cpu_t from_vector(const std::vector<data_t, allocator_t>&)
    {
        return cpu;
    }
    
    template <typename data_t, typename allocator_t>
    inline gpu_t from_vector(const device_vector<data_t, allocator_t>&)
    {
        return gpu;
    }
    
    //Might be best to do bounds checking here
    template <typename data_t, typename allocator_t>
    inline data_t inspect_at(const std::vector<data_t, allocator_t>& data, const std::size_t idx)
    {
        return data[idx];
    }
    
    template <typename data_t, typename allocator_t>
    inline data_t inspect_at(const device_vector<data_t, allocator_t>& data, const std::size_t idx)
    {
        data_t out;
#if (_sp_cuda)
        cudaMemcpy(&out, data.raw + idx, sizeof(data_t), cudaMemcpyDeviceToHost);
#endif
        return out;
    }
}