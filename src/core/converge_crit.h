#pragma once

#include <concepts>

namespace spade::utils
{
    template <std::floating_point float_t>
    struct converge_crit_t
    {
        float_t err_tol;
        int max_iters;
        
        float_t tolerance() const { return err_tol; }
        int     max_its()   const { return max_iters; }
    };
}