#pragma once

#include <concepts>

namespace spade::utils
{
    template <typename T, typename tag_t> concept has_tag = std::same_as<typename T::tag_type, tag_t>;
    
    template <typename tag_t, typename val_t = int>
    struct tagged_value_t
    {
        val_t value;
        using tag_type   = tag_t;
        using value_type = val_t;
        tagged_value_t(const val_t v) : value{v} {}
        operator val_t() const { return value; }
    };
    
    template <typename tag_t> struct tagger_t
    {
        template <typename val_t> tagged_value_t<tag_t, val_t> operator() (const val_t& v)
        {
            return tagged_value_t<tag_t, val_t>(v);
        }
    };
    
    struct tag_impl_t
    {
        template <typename tag_t> tagger_t<tag_t> operator[] (const tag_t& t) const
        {
            return tagger_t<tag_t>();
        }
    };
    static const tag_impl_t tag;
}