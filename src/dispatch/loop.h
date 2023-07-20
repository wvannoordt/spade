#pragma once

#include "dispatch/device_type.h"

namespace spade::dispatch
{
    template <typename range_t, typename loop_body_t>
    requires(device::is_cpu<typename range_t::device_t>)
    auto loop(const range_t& loop_range, const loop_body_t& body)
    {
        
    }
    
    template <typename range_t, typename loop_body_t>
    requires(device::is_gpu<typename range_t::device_t>)
    auto loop(const range_t& loop_range, const loop_body_t& body)
    {
        
    }
}