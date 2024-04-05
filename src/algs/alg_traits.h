#pragma once

#include "sym/sym.h"

namespace spade::algs
{
    template <typename T> concept is_trait_label = sym::is_symbol<T>;
    
    template <typename T> concept is_alg_trait = requires(T t)
    {
        { t.trait_label() } -> is_trait_label;
    };
    
    template <is_alg_trait... traits_t>
    struct trait_list_t
    {
        std::tuple<traits_t...> data;
    };
    
    struct not_found_t{};
    
    namespace detail
    {
        template <const int idx, is_trait_label key_t, is_alg_trait... traits_t>
        requires ((idx < 0) || idx >= sizeof...(traits_t))
        static auto search_get_trait(const trait_list_t<traits_t...>& list, const key_t& key)
        {
            return not_found_t();
        }
        
        template <const int idx, is_trait_label key_t, is_alg_trait... traits_t>
        requires ((idx >= 0) && idx < sizeof...(traits_t))
        static auto search_get_trait(const trait_list_t<traits_t...>& list, const key_t& key)
        {
            const auto& elem = std::get<idx>(list.data);
            const auto found_key = elem.trait_label();
            using found_key_t = typename utils::remove_all<decltype(found_key)>::type;
            if constexpr (std::same_as<found_key_t, key_t>)
            {
                return elem;
            }
            else
            {
                return search_get_trait<idx+1>(list, key);
            }
        }
    }
    
    template <is_trait_label key_t, is_alg_trait default_t, is_alg_trait... traits_t>
    static auto get_trait(const trait_list_t<traits_t...>& trait_list, const key_t& key, const default_t& def)
    {
        const auto def_label = def.trait_label();
        using def_label_t = typename utils::remove_all<decltype(def_label)>::type;
        static_assert(std::same_as<def_label_t, key_t>, "Default trait in get_trait must have the same label as provided in the search");
        const auto result = detail::search_get_trait<0>(trait_list, key);
        using result_type = typename utils::remove_all<decltype(result)>::type;
        constexpr bool return_default = std::same_as<result_type, not_found_t>;
        if constexpr (return_default)
        {
            return def;
        }
        else
        {
            return result;
        }
    }
    
    template <is_trait_label key_t, is_alg_trait... traits_t>
    static auto get_trait(const trait_list_t<traits_t...>& trait_list, const key_t& key)
    {
        const auto result = detail::search_get_trait<0>(trait_list, key);
        using result_type = typename utils::remove_all<decltype(result)>::type;
        static_assert(!std::same_as<result_type, not_found_t>, "Could not find the required trait in the provided list");
        return result;
    }
    
    template <is_alg_trait... traits_t>
    static auto make_traits(const traits_t&... traits)
    {
        return trait_list_t<traits_t...>{std::tuple<traits_t...>{traits...}};
    }
}