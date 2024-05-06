#pragma once

#include <concepts>

#include "grid/grid_index_types.h"
#include "core/utils.h"
#include "grid/grid.h"

namespace spade::omni
{
    namespace info
    {
        template <typename T> concept is_omni_info = requires(const T& t) {T::info_tag;};
        template <typename derived_t> struct info_base
        {
            constexpr static bool info_tag = true;
            
            std::string get_name() const
            {
                return "[base]";
            }
            
            std::string name() const
            {
                return static_cast<derived_t*>(this)->get_name();
            }
        };

        static struct undirected_tag_t{} undirected;

        consteval static bool is_direction(const int& i)              {return true;}
        consteval static bool is_direction(const undirected_tag_t& i) {return false;}
    }

    template <typename... infos_t>
    requires (info::is_omni_info<infos_t> && ...)
    struct info_list_t
    {
        constexpr static int num_infos() {return sizeof...(infos_t);}
        template <typename info_t> constexpr static bool contains = (std::same_as<info_t, infos_t> || ...);
        template <const int idx>   using info_elem = typename utils::get_pack_type<idx, infos_t...>::type;
        template <typename add_t>  using expand_list = info_list_t<infos_t..., add_t>;
    };
}