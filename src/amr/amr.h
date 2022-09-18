#pragma once

#include "amr/amr_blocks.h"

namespace spade::amr
{
    template <typename blocks_t, typename anode_t>
    void refine(
        blocks_t& blocks,
        anode_t& anode,
        const typename anode_t::amr_refine_t& refine_mode = typename anode_t::amr_refine_t(true, true, true))
    {
        anode.refine_node_recurse(refine_mode);
        blocks.enumerate();
    }
}