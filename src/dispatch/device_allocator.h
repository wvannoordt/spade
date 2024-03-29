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
    struct device_allocator_t
    {
        int count = 0;
        using value_type = data_t;
        device_allocator_t() = default;
        
        template<typename rhs_t>
        constexpr device_allocator_t (const device_allocator_t <rhs_t>&) noexcept {}
        
        constexpr static std::size_t max_alloc_size = std::numeric_limits<std::size_t>::max() / sizeof(data_t);
        [[nodiscard]] data_t* allocate(std::size_t n)
        {
            if (n > max_alloc_size) throw mem_exception("attempted to allocate too much device memory (" + std::to_string(n*sizeof(data_t)) + " bytes)");
            data_t* output = nullptr;
            if (n==0) return output;
            std::string err_string = "attempted to allocate device memory without device support";
#if(_sp_cuda)
            
            auto er_code = cudaMalloc(&output, n*sizeof(data_t));
            if (er_code == cudaSuccess && output != nullptr)
            {
                ++count;
                return output;
            }
            err_string = std::string(cudaGetErrorString(er_code));
#endif
            throw mem_exception("unable to allocate on device: " + err_string);
        }
        
        void deallocate(data_t* p, std::size_t n) noexcept
        {
            if (p == nullptr) return;
#if(_sp_cuda)
            --count;
            auto er_code = cudaFree(p);
            
            //Done to avoid annoying compiler warning.
            p = nullptr;
            if (er_code != cudaSuccess)
            {
                throw mem_exception("could not deallocate device memory: " + std::string(cudaGetErrorString(er_code)));
            }
#endif
        }
    };

#if(_sp_cuda)
    template <typename data_t>
    struct pinned_allocator_t
    {
        int count = 0;
        using value_type = data_t;
        pinned_allocator_t() = default;
        
        template<typename rhs_t>
        constexpr pinned_allocator_t (const pinned_allocator_t <rhs_t>&) noexcept {}
        
        constexpr static std::size_t max_alloc_size = std::numeric_limits<std::size_t>::max() / sizeof(data_t);
        [[nodiscard]] data_t* allocate(std::size_t n)
        {
            if (n > max_alloc_size) throw mem_exception("attempted to allocate too much memory (" + std::to_string(n*sizeof(data_t)) + " bytes)");
            data_t* output = nullptr;
            if (n==0) return output;
            std::string err_string = "attempted to allocate host-pinned memory without device support";
            auto er_code = cudaMallocHost(&output, n*sizeof(data_t), cudaHostAllocDefault);
            if (er_code == cudaSuccess && output != nullptr)
            {
                ++count;
                return output;
            }
            err_string = std::string(cudaGetErrorString(er_code));
            throw mem_exception("unable to allocate pinned memory: " + err_string);
        }
        
        void deallocate(data_t* p, std::size_t n) noexcept
        {
            if (p == nullptr) return;
            --count;
            auto er_code = cudaFreeHost(p);
            
            //Done to avoid annoying compiler warning.
            p = nullptr;
            if (er_code != cudaSuccess)
            {
                throw mem_exception("could not deallocate pinned memory: " + std::string(cudaGetErrorString(er_code)));
            }
        }
    };
#else
    template <typename data_t> using pinned_allocator_t = std::allocator<data_t>;
#endif
}