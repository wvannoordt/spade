#pragma once

#include "core/config.h"
#include "core/ctrs.h"
#include "core/bounding_box.h"

namespace spade::algs
{
    template <
        typename index_value_t,
        const std::size_t index_rank,
        typename fcn_t,
        const int index_val,
        typename index_t
        >
    requires (index_val >= 0)
    static void exec_loop(
        const spade::bound_box_t<index_value_t, index_rank>& bds,
        index_t& index,
        const fcn_t& fcn)
    {
        auto imin = bds.min(index_val);
        auto imax = bds.max(index_val);
        for (index_value_t i = imin; i < imax; ++i)
        {
            index[index_val] = i;
            exec_loop<index_value_t, index_rank, fcn_t, index_val-1, index_t>(bds, index, fcn);
        }
    }

    template <
        typename index_value_t,
        const std::size_t index_rank,
        typename fcn_t,
        const int index_val,
        typename index_t
        >
    requires (index_val < 0)
    static void exec_loop(
        const spade::bound_box_t<index_value_t, index_rank>& bds,
        const index_t& index,
        const fcn_t& fcn)
    {
        fcn(index);
    }

    template <
        typename index_value_t,
        const std::size_t index_rank,
        typename fcn_t,
        typename index_t = spade::ctrs::array<index_value_t, index_rank>
        >
    static void md_loop(
        const spade::bound_box_t<index_value_t, index_rank>& bds,
        const fcn_t& fcn)
    {
        index_t index;
        exec_loop<index_value_t, index_rank, fcn_t, index_rank-1, index_t>(bds, index, fcn);
    }
}