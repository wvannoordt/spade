#pragma once

#include "dispatch/ranges/grid_idx_range.h"
#include "core/bounding_box.h"

namespace spade::dispatch
{
    template <typename arr_t>
    auto support_of(const arr_t& arr, const grid::exchange_inclusion_e& i_exchg)
    {
        using index_t    = typename arr_t::index_type;
        static_assert(index_t::centering_type() == grid::cell_centered, "support_of has no convention for face, node, or edge cenetering yet!");
        const auto& grid = arr.get_grid();
        const auto dev   = arr.device();
        index_t low;
        
        int iexch = (i_exchg == grid::include_exchanges)?1:0;
        
        low.i()  = 0 - iexch*arr.get_num_exchange(0);
        low.j()  = 0 - iexch*arr.get_num_exchange(1);
        low.k()  = 0 - iexch*arr.get_num_exchange(2);
        low.lb() = 0;
        
        index_t high;
        high.i()  = grid.get_num_cells(0) + iexch*arr.get_num_exchange(0);
        high.j()  = grid.get_num_cells(1) + iexch*arr.get_num_exchange(1);
        high.k()  = grid.get_num_cells(2) + iexch*arr.get_num_exchange(2);
        high.lb() = grid.get_num_local_blocks();
        
        if constexpr (grid::is_directed_index<index_t>)
        {
            low.dir()  = 0;
            high.dir() = grid.dim();
        }
        
        return ranges::grid_idx_range_t(low, high, dev);
    }
    
    template <typename arr_t>
    auto support(const arr_t& arr, const grid::exchange_inclusion_e& i_exchg)
    {
        using index_t    = typename arr_t::index_type;
        static_assert(index_t::centering_type() == grid::cell_centered, "support_of has no convention for face, node, or edge cenetering yet!");
        const auto& grid = arr.get_grid();
        const auto dev   = arr.device();
        index_t low;
        
        int iexch = (i_exchg == grid::include_exchanges)?1:0;
        
        low.i()  = 0 - iexch*arr.get_num_exchange(0);
        low.j()  = 0 - iexch*arr.get_num_exchange(1);
        low.k()  = 0 - iexch*arr.get_num_exchange(2);
        low.lb() = 0;
        
        index_t high;
        high.i()  = grid.get_num_cells(0) + iexch*arr.get_num_exchange(0);
        high.j()  = grid.get_num_cells(1) + iexch*arr.get_num_exchange(1);
        high.k()  = grid.get_num_cells(2) + iexch*arr.get_num_exchange(2);
        high.lb() = grid.get_num_local_blocks();
        
        if constexpr (grid::is_directed_index<index_t>)
        {
            low.dir()  = 0;
            high.dir() = grid.dim();
        }
        
        auto bnds = utils::make_bounds(low.i(), high.i(), low.j(), high.j(), low.k(), high.k(), low.lb(), high.lb());
        ranges::basic_range_t<typename index_t::value_type, index_t::size(), index_t> rg{bnds};
        return rg;
    }
}