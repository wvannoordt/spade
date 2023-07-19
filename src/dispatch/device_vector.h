#pragma once

#include "dispatch/device_allocator.h"

namespace spade::device
{
    //todo: improve this implementation
    template <typename data_t>
    struct vector
    {
        using allocator_t = cuda_allocator_t<data_t>;
        allocator_t allocator;
        data_t* raw = nullptr;
        std::size_t c_size = 0;
        
        vector(){}
        vector(const std::size_t n) 
        {
            this->resize(n);
        }
        
        void resize(const std::size_t& n)
        {
            if (raw != nullptr) allocator.deallocate(raw, c_size);
            c_size = n;
            raw = allocator.allocate(n);
        }
        
        ~vector()
        {
            if (raw != nullptr) allocator.deallocate(raw, c_size);
        }
        
        _sp_hybrid data_t& operator[] (const std::size_t& idx)             { return raw[idx]; }
        _sp_hybrid const data_t& operator[] (const std::size_t& idx) const { return raw[idx]; }
    };
}