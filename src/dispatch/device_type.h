#pragma once

#include <type_traits>
#include <concepts>
#include "core/cuda_incl.h"
namespace spade::device
{
    struct gpu_t {};
    static gpu_t gpu;
    
    struct cpu_t {};
    static cpu_t cpu;
    
    using best_type = std::conditional<_sp_cuda, gpu_t, cpu_t>::type;
    static best_type best;
    
    template <typename T> concept is_cpu    = std::same_as<T, cpu_t>;
    template <typename T> concept is_gpu    = std::same_as<T, gpu_t>;
    
    template <typename T> concept is_device = is_cpu<T> || is_gpu<T>;
    
    template <typename T> concept has_device = requires(const T& t)
    {
        t.device();
    };
}
