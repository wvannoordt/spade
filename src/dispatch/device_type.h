#pragma once

#include <concepts>

namespace spade::device
{
    struct gpu_t {} gpu;
    struct cpu_t {} cpu;
    
  template <typename T> concept is_cpu = std::same_as<T, cpu_t>;
  template <typename T> concept is_gpu = std::same_as<T, gpu_t>;
}
