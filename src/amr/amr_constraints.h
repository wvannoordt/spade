#pragma once

#include "amr_node.h"

namespace spade::amr
{
    namespace constraints
    {
        const static auto factor2 = [](const auto& node, const auto& neigh)
        {
            using out_t = typename utils::remove_all<decltype(node)>::type::amr_refine_t;
            const auto& neigh_node = neigh.endpoint.get();
            out_t out = false;
            for (int i = 0; i < out.size(); ++i)
            {
                if (node.level[i] > neigh_node.level[i] + 1) out[i] = true;
            }
            return out;
        };
    }
}