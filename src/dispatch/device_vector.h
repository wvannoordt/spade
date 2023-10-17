#pragma once

#include "dispatch/device_allocator.h"

namespace spade::device
{
    //todo: improve this implementation
    template <typename data_t>
    struct device_vector
    {
        using value_type = data_t;
        using allocator_t = device_allocator_t<data_t>;
        allocator_t allocator;
        data_t* raw = nullptr;
        std::size_t c_size = 0;
        
        device_vector(){}
        device_vector(const std::size_t n) 
        {
            this->resize(n);
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
        
        std::size_t size() const {return c_size;}
        
        void resize(const std::size_t& n)
        {
#if (_sp_cuda)
            if (raw != nullptr) allocator.deallocate(raw, c_size);
            c_size = n;
            raw = allocator.allocate(n);
#endif
        }
        
        ~device_vector()
        {
            if (raw != nullptr) allocator.deallocate(raw, c_size);
        }
        
        _sp_hybrid data_t& operator[]       (const std::size_t& idx)       { return raw[idx]; }
        _sp_hybrid const data_t& operator[] (const std::size_t& idx) const { return raw[idx]; }
    };
}