#pragma once

#include "core/cuda_incl.h"

namespace spade::device
{
    struct mem_exception : public std::exception
    {
        std::string message;
        mem_exception(const std::string& message_in) : message{message_in} {}
        const char* what() const throw()
        {
            return message.c_str();
        }
    };
    
    template <typename data_t>
    struct cuda_allocator_t
    {
        int count = 0;
        using value_type = data_t;
        cuda_allocator_t() = default;
        
        template<typename rhs_t>
        constexpr cuda_allocator_t (const cuda_allocator_t <rhs_t>&) noexcept {}
        
        constexpr static std::size_t max_alloc_size = std::numeric_limits<std::size_t>::max() / sizeof(data_t);
        [[nodiscard]] data_t* allocate(std::size_t n)
        {
            if (n > max_alloc_size) throw mem_exception("attempted to allocate too much device memory (" + std::to_string(n*sizeof(data_t)) + " bytes)");
            data_t* output = nullptr;
            void** alias = (void**)&output;
            auto er_code = cudaMalloc(alias, n*sizeof(data_t));
            if (er_code == cudaSuccess && output != nullptr)
            {
                ++count;
                return output;
            }
            throw mem_exception("unable to allocate on device: " + std::string(cudaGetErrorString(er_code)));
        }
        
        void deallocate(data_t* p, std::size_t n) noexcept
        {
            --count;
            cudaFree((void*)p);
        }
    };
}