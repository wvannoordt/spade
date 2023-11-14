#pragma once

#include "sym/symbol.h"

namespace spade::ode
{
    template <sym::is_symbol sym_t, typename func_t>
    struct expression_t
    {
        using symbol_type = sym_t;
        const func_t func;
        expression_t(const sym_t&, const func_t& f_in)
        : func{f_in} {}
    };
}