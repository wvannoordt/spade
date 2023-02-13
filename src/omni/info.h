#pragma once

#include "core/grid_index_types.h"

namespace spade::omni
{
    namespace info
    {
        template <typename derived_t> struct info_base
        {
            constexpr static bool info_tag = true;
        };
        struct index    : public info_base<index    > {};
        struct normal   : public info_base<normal   > {}; // only supported at locations centered at faces, expressed in index space
        struct coord    : public info_base<coord    > {};
        struct metric   : public info_base<metric   > {};
        struct jacobian : public info_base<jacobian > {};
        struct value    : public info_base<value    > {}; // interpolation to nodes?
        struct gradient : public info_base<gradient > {}; // floating-point types only, and need to be able to specify the order-of-accuracy later
    }
    template <typename... infos_t>
    requires (infos_t::info_tag && ...)
    struct info_list_t {};
}