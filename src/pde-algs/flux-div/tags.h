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
    
    static struct tbasic_t : public fdiv_alg_base_t {} basic;
    static struct tbfoct_t : public fdiv_alg_base_t {} bfoct;
    static tldbal_t<false>    ldbalnp;
    static tldbal_t<true>     ldbal;
}