#pragma once

#include "core/ctrs.h"

namespace spade::algs
{
    template <typename dest_t, typename src_t> auto&       unsafe_cast(src_t& src)       {return *((dest_t*)(&src));}
    template <typename dest_t, typename src_t> const auto& unsafe_cast(const src_t& src) {return *((dest_t*)(&src));}
}