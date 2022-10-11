#pragma once

#include "amr/amr_blocks.h"

namespace spade::amr
{
    template <typename blocks_t>
    void refine(
        blocks_t& blocks,
        const std::size_t& index,
        const typename blocks_t::node_type::amr_refine_t& refine_mode = typename blocks_t::node_type::amr_refine_t(true, true, true))
    {
        typename blocks_t::node_type& anode = *(blocks.all_nodes[index]);
        anode.refine_node_recurse(refine_mode);
        blocks.enumerate();
    }
}