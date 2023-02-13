#pragma once

#include "core/grid_index_types.h"

namespace spade::omni
{
    template <const int di, const int dj, const int dk> struct offset_t
    {
        constexpr static bool is_1d = ((dj==0) && (di==0));
    };
}