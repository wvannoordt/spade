#pragma once

#include "grid/grid.h"
#include "dispatch/support_of.h"
#include "dispatch/device_type.h"
#include "dispatch/execute.h"

namespace spade::algs
{
    template <grid::multiblock_array array_t, class callable_t>
    requires (std::invocable<callable_t, typename array_t::index_type>)
    void for_each(array_t& arr, callable_t& func, const grid::exchange_inclusion_e& exchange_policy = grid::exclude_exchanges)
    {
        using index_type        = typename array_t::index_type;
        using val_t             = typename array_t::alias_type;
        constexpr int tile_size = 4;
        const auto tile_range   = dispatch::ranges::make_range(0, tile_size, 0, tile_size, 0, tile_size);
        dispatch::kernel_threads_t kpool(tile_range, arr.device());
        using threads_type = decltype(kpool);
        
        auto nx        = arr.get_grid().get_num_cells();
        auto ng        = arr.get_num_exchange();
        auto nx_extent = nx;
        auto ntiles    = nx_extent;
        if (exchange_policy == grid::include_exchanges)
        {
            for (int d = 0; d < arr.get_grid().dim(); ++d)
            {
                nx_extent[d] += ng[d]*2;
            }
        }
        
        for (int d = 0; d < nx.size(); ++d)
        {
            if ((nx[d] % tile_size) != 0)
            {
                throw except::sp_exception("transform_inplace requires a block size that is a multiple of 4");
            }
            if (ng[d] != 2)
            {
                throw except::sp_exception("transform_inplace requires all exchange cells are size 2");
            }
            ntiles[d]    = utils::i_div_up(nx_extent[d], tile_size);
        }
        
        const auto outer_range = dispatch::ranges::make_range(0, ntiles[0], 0, ntiles[1]*ntiles[2], 0, int(arr.get_grid().get_num_local_blocks()));
        int iexch = int(exchange_policy == grid::include_exchanges);
        auto loop_load = [=] _sp_hybrid (const ctrs::array<int, 3>& outer_raw, const threads_type& threads) mutable
        {
            int tile_id_1d = outer_raw[1];
            ctrs::array<int, 3> tile_id;
            tile_id[0] = outer_raw[0];
            // tile_id[0]  = tile_id_1d % ntiles[0];
            // tile_id_1d -= tile_id[0];
            // tile_id_1d /= ntiles[0];
            tile_id[1]  = tile_id_1d % ntiles[1];
            tile_id_1d -= tile_id[1];
            tile_id_1d /= ntiles[1];
            tile_id[2]  = tile_id_1d;
            threads.exec([&](const ctrs::array<int, 3>& inner_raw)
            {
                int i  = tile_id[0]*tile_size + inner_raw[0] - iexch*ng[0];
                int j  = tile_id[1]*tile_size + inner_raw[1] - iexch*ng[1];
                int k  = tile_id[2]*tile_size + inner_raw[2] - iexch*ng[2];
                int lb = outer_raw[2];
                index_type index(i, j, k, lb);
                func(index);
            });
        };
        dispatch::execute(outer_range, loop_load, kpool);
    }
}