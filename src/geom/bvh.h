#pragma once

#include <string>

namespace spade::geom
{
    template <const int dim, typename handle_t, typename float_t, template <typename> typename container_tt>
    struct bvh_t
    {
        using idx_t = std::size_t;
        static constexpr idx_t no_val = std::string::npos;
        
        //idea; we will organize the stored handles by a depth-first traversal
        container_tt<handle_t> data;
        container_tt<idx_t> children;
    };
}