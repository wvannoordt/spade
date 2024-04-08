#pragma once

#include <vector>
#include <algorithm>
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
    
    inline int device_count()
    {
        int out = 0;
#if (_sp_cuda)
        cudaGetDeviceCount(&out);
#endif
        return out;
    }
    
    inline void set_device(const int dev_id)
    {
#if (_sp_cuda)
        auto ercode = cudaSetDevice(dev_id);
        if (ercode != cudaSuccess)
        {
            throw except::sp_exception("Attempted to set an invalid device");
        }
#endif
    }
    
    inline bool enable_p2p(int thread_id, std::vector<int> devices)
    {
        // std::sort(devices.begin(), devices.end());
        // if (std::unique(devices.begin(), devices.end()) != devices.end())
        // {
        //     throw except::sp_exception("Duplicate device in launch configuration!");
        // }
#if (_sp_cuda)
        for (auto d: devices)
        {
            if (d != devices[thread_id])
            {
                int can_access;
                cudaDeviceCanAccessPeer(&can_access, devices[thread_id], d);
                if (can_access != 1)
                {
                    return false;
                }
                else
                {
                    auto code = cudaDeviceEnablePeerAccess(d, 0);
                    if (code == cudaErrorInvalidDevice)
                    {
                        throw except::sp_exception("Failed to enable P2P access on device " + std::to_string(d) + ": invalid device");
                    }
                    if (code == cudaErrorInvalidValue)
                    {
                        throw except::sp_exception("Failed to enable P2P access on device " + std::to_string(d) + ": invalid value");
                    }
                }
            }
        }
#endif
        return true;
    }
}
