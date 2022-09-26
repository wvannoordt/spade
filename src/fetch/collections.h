#pragma once

#include <tuple>
#include <type_traits>
#include <concepts>
#include <typeinfo>

#include "core/grid.h"
#include "core/coord_system.h"

namespace spade::fetch
{
    template <typename... infos_t> struct cell_info
    {
        const static std::size_t num_params = sizeof...(infos_t);
        std::tuple<infos_t...> elements;
        template <const std::size_t idx> const auto& item(void) const {return std::get<idx>(elements).data;}
        template <const std::size_t idx> auto& item(void) {return std::get<idx>(elements).data;}
    };

    template <typename... infos_t>
    static std::ostream & operator<<(std::ostream & os, const cell_info<infos_t...>& ftch)
    {
        static_for<0,ftch.num_params>([&](const auto& i) -> void
        {
            os << "Data Item " << i.value << ":\n" << std::get<i.value>(ftch.elements) << "\n";
        });
        return os;
    }

    template <typename... infos_t> struct face_info
    {
        const static std::size_t num_params = sizeof...(infos_t);
        std::tuple<infos_t...> elements;
        template <const std::size_t idx> const auto& item(void) const {return std::get<idx>(elements).data;}
        template <const std::size_t idx> auto& item(void) {return std::get<idx>(elements).data;}
    };

    template <typename... infos_t>
    static std::ostream & operator<<(std::ostream & os, const face_info<infos_t...>& ftch)
    {
        static_for<0,ftch.num_params>([&](const auto& i) -> void
        {
            os << "Data Item " << i.value << ":\n" << std::get<i.value>(ftch.elements) << "\n";
        });
        return os;
    }
}