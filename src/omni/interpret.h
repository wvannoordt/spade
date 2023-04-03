#pragma once

#include "omni/info.h"
#include "omni/stencil.h"
#include "omni/info_union.h"
#include "omni/stencil_union.h"

namespace spade::omni
{
    //Note that input_t could be an alias itself!
    template <typename interpret_t, typename input_t>
    struct stencil_alias_t
    {
        using stencil_type        = interpret_t; // used for any further (deeper) aliases
        using native_stencil_type = typename input_t::stencil_type;
        const input_t& native;
        stencil_alias_t(const input_t& native_in) : native{native_in}{}
    };

    template <typename interpret_t, typename input_t>
    auto interpret_stencil(const input_t& input)
    {
        return stencil_alias_t<interpret_t, input_t>(input);
    }
}