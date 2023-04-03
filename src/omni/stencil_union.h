#pragma once

#include "omni/stencil.h"

namespace spade::omni
{
    namespace detail
    {
        template <typename stencil1_t, typename stencil2_t> struct stencil_union_t
        {
            using type = int;
        };
    }

    template <typename stencil1_t, typename stencil2_t> using stencil_union = detail::stencil_union_t<stencil1_t, stencil2_t>::type;
}