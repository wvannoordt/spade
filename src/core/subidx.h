#pragma once

#include "core/grid_index_types.h"

namespace spade::subidx
{
    template <typename index_t> auto i(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::i_subindex>::value];
    }
    
    template <typename index_t> auto j(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::j_subindex>::value];
    }
    
    template <typename index_t> auto k(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::k_subindex>::value];
    }
    
    template <typename index_t> auto lb(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::lb_subindex>::value];
    }
    
    template <typename index_t> auto dir(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering, grid::dir_subindex>::value];
    }
}