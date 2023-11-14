#pragma once

#include "core/cuda_incl.h"
#include "core/ctrs.h"

namespace spade::algs
{
    // NOTE!!!!!!!!!!!!!!
    // We need to re-write this using bit_cast!!!!!!!!!!!
    template <typename dest_t, typename src_t> _sp_hybrid auto&       unsafe_cast(src_t& src)       {return *((dest_t*)(&src));}
    template <typename dest_t, typename src_t> _sp_hybrid const auto& unsafe_cast(const src_t& src) {return *((dest_t*)(&src));}
}