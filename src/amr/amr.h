#pragma once

#include "amr/amr_blocks.h"
#include "amr/amr_constraints.h"
namespace spade::amr
{
    template <typename blocks_t>
    void refine(
        blocks_t& blocks,
        const std::size_t& index,
        const typename blocks_t::node_type::amr_refine_t& refine_mode = typename blocks_t::node_type::amr_refine_t(true, true, true))
    {
        
    }
}