#pragma once

#include "core/ctrs.h"

namespace spade::dispatch::shmem
{
    struct empty_t
    {
        constexpr static std::size_t bytes() { return 0; }
    };
}