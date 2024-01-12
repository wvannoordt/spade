#pragma once

#include "core/ctrs.h"

namespace spade::dispatch::shmem
{
    struct empty_t
    {
        constexpr static std::size_t bytes() { return 0; }
    };
    
    template<typename data_t, typename... datas_t>
    struct shmem_t
    {
        char* base;
    };
    
    template <typename... datas_t>
    auto make_shmem(const datas_t&... datas)
    {
        
    }
}