#pragma once

#include <concepts>

#include "core/bounding_box.h"

namespace spade::dispatch::ranges
{
    template <std::integral idx_t, const std::size_t ar_size, typename index_t = ctrs::array<idx_t, ar_size>>
    struct basic_range_t
    {
        bound_box_t<idx_t, ar_size> bounds;
        
        using index_type = index_t;
        
        _sp_hybrid bool is_empty() const { return bounds.volume() == 0; }
    };
    
    template <std::integral idx_t, const std::size_t rank>
    inline auto make_range(const bound_box_t<idx_t, rank>& bounds)
    {
        return basic_range_t<idx_t, rank>{bounds};
    }
    
    template <std::integral... idxs_t>
    inline auto make_range(const idxs_t&... vs)
    {
        auto bnd = utils::make_bounds(vs...);
        return make_range(bnd);
    }
    
    template <std::integral idx_t, typename index_type>
    _sp_hybrid index_type compute_index(
        const basic_range_t<idx_t, std::size_t(1), index_type>& i_range,
        const ranges::outer_config_t& g_dim,
        const ranges::outer_config_t& b_idx)
    {
        return b_idx.data.x;
    }
    
    template <std::integral idx_t, typename index_type>
    _sp_hybrid index_type compute_index(
        const basic_range_t<idx_t, std::size_t(1), index_type>& i_range,
        const ranges::inner_config_t& b_dim,
        const ranges::inner_config_t& t_idx)
    {
        return t_idx.data.x;
    }
    
    template <std::integral idx_t, typename index_type>
    _sp_hybrid index_type compute_index(
        const basic_range_t<idx_t, std::size_t(1), index_type>& i_range,
        const ranges::outer_config_t& g_dim,
        const ranges::outer_config_t& b_idx,
        const ranges::inner_config_t& b_dim,
        const ranges::inner_config_t& t_idx)
    {
        return b_idx.data.x*b_dim.data.x + t_idx.data.x;
    }
    
    template <std::integral idx_t, typename index_type>
    _sp_hybrid index_type compute_index(
        const basic_range_t<idx_t, std::size_t(3), index_type>& i_range,
        const ranges::outer_config_t& g_dim,
        const ranges::outer_config_t& b_idx)
    {
        return {b_idx.data.x-i_range.bounds.min(0), b_idx.data.y-i_range.bounds.min(1), b_idx.data.z-i_range.bounds.min(2)};
    }
}