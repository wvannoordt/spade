#pragma once

#include "sym/sym.h"

namespace spade::pde_algs
{
    struct fdiv_alg_base_t
    {
        constexpr static auto trait_label() { using namespace sym::literals; return "flux_div"_sym; }
    };
    
    template <const bool v_use_parity>
    struct tldbal_t : public fdiv_alg_base_t
    {
        constexpr static bool use_parity = v_use_parity;
    };
    
    template <const bool v_use_simple_buffering>
    struct tfldbc_t : public fdiv_alg_base_t
    {
        constexpr static bool use_simple_buffering = v_use_simple_buffering;
    };
    
    static struct tbasic_t : public fdiv_alg_base_t {} basic;
    static struct tbfoct_t : public fdiv_alg_base_t {} bfoct;
    static struct tfused_t : public fdiv_alg_base_t {} fused;
    static tldbal_t<false>                             ldbalnp;
    static tldbal_t<true>                              ldbal;
    static tfldbc_t<false> fldbc;
    static tfldbc_t<true>  fldbcsb;
    
}