#pragma once

#include "core/grid_index_types.h"
#include "core/utils.h"
#include "core/grid.h"

namespace spade::omni
{
    namespace info
    {
        template <typename T> concept is_omni_info = requires(const T& t) {T::info_tag;};
        template <typename derived_t> struct info_base
        {
            constexpr static bool info_tag = true;
        };
    }

    template <typename... infos_t>
    requires (info::is_omni_info<infos_t> && ...)
    struct info_list_t
    {
        constexpr static int num_infos() {return sizeof...(infos_t);}

        template <const int idx>
        using info_elem = typename utils::get_pack_type<idx, infos_t...>::type;
    };
}