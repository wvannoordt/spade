#pragma once

#include <concepts>

#include "core/ctrs.h"
#include "core/bounding_box.h"

namespace spade::algs
{
    template <
        typename index_t,
        const int index_rank,
        typename fcn_t,
        const int index_val
        >
    requires (index_val >= 0)
    static void exec_loop(
        const spade::bound_box_t<index_t, index_rank>& bds,
        spade::ctrs::array<index_t, index_rank>& index,
        const fcn_t& fcn)
    {
        auto imin = bds.min(index_val);
        auto imax = bds.max(index_val);
        for (index_t i = imin; i < imax; ++i)
        {
            index[index_val] = i;
            exec_loop<index_t, index_rank, fcn_t, index_val-1>(bds, index, fcn);
        }
    }

    template <
        typename index_t,
        const int index_rank,
        typename fcn_t,
        const int index_val
        >
    requires (index_val < 0)
    static void exec_loop(
        const spade::bound_box_t<index_t, index_rank>& bds,
        const spade::ctrs::array<index_t, index_rank>& index,
        const fcn_t& fcn)
    {
        fcn(index);
    }

    template <
        typename index_t,
        const std::size_t index_rank,
        typename fcn_t
        >
    static void multi_loop(
        const spade::bound_box_t<index_t, index_rank> bds,
        const fcn_t& fcn)
    {
        spade::ctrs::array<index_t, index_rank> index;
        exec_loop<index_t, index_rank, fcn_t, index_rank-1>(bds, index, fcn);
    }
}