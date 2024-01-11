#pragma once

#include <concepts>

#include "core/bounding_box.h"

namespace spade::dispatch::ranges
{
    template <std::integral idx_t, const std::size_t ar_size>
    struct basic_range_t
    {
        bound_box_t<idx_t, ar_size> bounds;
        
        using index_type = ctrs::array<idx_t, ar_size>;
        
        bool is_empty() const { return bounds.volume() == 0; }
    };
    
    template <std::integral idx_t, const std::size_t rank>
    inline auto make_range(const bound_box_t<idx_t, rank>& bounds)
    {
        return basic_range_t<idx_t, rank>{bounds};
    }
}