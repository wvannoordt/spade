#pragma once

#include "core/attribs.h"
#include "core/grid_index_types.h"

namespace spade::subidx
{
    template <typename index_t> static _finline_ auto i(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::i_subindex>::value];
    }
    
    template <typename index_t> static _finline_ auto j(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::j_subindex>::value];
    }
    
    template <typename index_t> static _finline_ auto k(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::k_subindex>::value];
    }
    
    template <typename index_t> static _finline_ auto lb(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering,  grid::lb_subindex>::value];
    }
    
    template <typename index_t> static _finline_ auto dir(const index_t& idx)
    {
        const auto centering = index_t::value_type::array_centering();
        return idx[grid::get_subidx<centering, grid::dir_subindex>::value];
    }
}