#pragma once

namespace spade::pde_algs
{
    
    template <const bool v_use_parity>
    struct tldbal_t
    {
        constexpr static bool use_parity = v_use_parity;
    };
    
    static struct tbasic_t {} basic;
    static struct tbfoct_t {} bfoct;
    static tldbal_t<false>    ldbalnp;
    static tldbal_t<true>     ldbal;
}