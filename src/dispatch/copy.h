#pragma once

#include <concepts>
#include "core/cuda_incl.h"

namespace spade::device
{
    // bad
    template <typename T> concept is_buffer = requires(const T t)
    {
        t.size();
        typename T::value_type;
        &t[0];
    };
    
    template <is_buffer send_t, is_buffer recv_t>
    requires (std::same_as<typename send_t::value_type, typename recv_t::value_type>)
    static void copy(const send_t& source, recv_t& dest)
    {
        const std::size_t send_size = sizeof(typename send_t::value_type)*source.size();
#if (_sp_cuda)
        auto er_code = cudaMemcpy(&dest[0], &source[0], send_size, cudaMemcpyDeviceToDevice);
        if (er_code != cudaSuccess)
        {
            std::string msg = "device copy failed: ";
            msg += cudaGetErrorString(er_code);
            throw except::sp_exception(msg);
        }
#else
        throw except::sp_exception("attempted device copy without gpu support");
#endif
    }
    
    template <is_buffer send_t, is_buffer recv_t>
    requires (std::same_as<typename send_t::value_type, typename recv_t::value_type>)
    static void peer_copy(const send_t& source, int sender, recv_t& dest, int receiver)
    {
        const std::size_t send_size = sizeof(typename send_t::value_type)*source.size();
#if (_sp_cuda)
        auto er_code = cudaMemcpyPeer(&dest[0], receiver, &source[0], sender, send_size);
        if (er_code != cudaSuccess)
        {
            std::string msg = "peer device copy failed: ";
            msg += cudaGetErrorString(er_code);
            throw except::sp_exception(msg);
        }
#else
        throw except::sp_exception("attempted device peer copy without gpu support");
#endif
    }
}