#pragma once

#include "grid/grid.h"
#include "dispatch/support_of.h"
#include "dispatch/device_type.h"
#include "dispatch/execute.h"
#include "omni/omni.h"
#include "core/invoke.h"

namespace spade::algs
{
    
    // template <grid::multiblock_array array_t, class callable_t>
    // void transform_inplace(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy = grid::exclude_exchanges)
    // {
    //     const auto ctr       = array_t::centering_type();
    //     auto kernel          = omni::to_omni<ctr>(func, arr);
    //     auto var_range       = dispatch::support_of(arr, exchange_policy);
    //     using index_type     = typename decltype(var_range)::index_type;
    //     auto d_image         = arr.image();
    //     const auto g_image   = arr.get_grid().image(arr.device());
        
    //     using val_t = typename array_t::alias_type;
    //     auto loop_load = [=] _sp_hybrid (const index_type& index) mutable
    //     {
    //         const val_t data = invoke_at(g_image, d_image, index, kernel);
    //         d_image.set_elem(index, data);
    //     };
    //     dispatch::execute(var_range, loop_load);
    // }
    
    // This seems to be marginally better than the simple one above
    template <grid::multiblock_array array_t, class callable_t>
    void transform_inplace(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy = grid::exclude_exchanges)
    {
        const auto ctr          = array_t::centering_type();
        auto kernel             = omni::to_omni<ctr>(func, arr);
        auto d_image            = arr.image();
        const auto g_image      = arr.get_grid().image(arr.device());
        using val_t             = typename array_t::alias_type;
        constexpr int tile_size = 4;
        using index_type = typename array_t::index_type;
        
        constexpr int dim = array_t::grid_type::dim();
        constexpr bool is_3d = dim == 3;
        
        ctrs::array<int, 3> ts{tile_size, tile_size, tile_size};
        if constexpr (!is_3d)
        {
            ts = {2*tile_size, 2*tile_size, 1};
        }
        const auto tile_range   = dispatch::ranges::make_range(0, ts[0], 0, ts[1], 0, ts[2]);
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
        
        ntiles = 1;
        for (int d = 0; d < arr.get_grid().dim(); ++d)
        {
            if ((nx[d] % tile_size) != 0)
            {
                throw except::sp_exception("transform_inplace requires a block size that is a multiple of 4");
            }
            ntiles[d] = utils::i_div_up(nx_extent[d], ts[d]);
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
                int i  = tile_id[0]*ts[0] + inner_raw[0] - iexch*ng[0];
                int j  = tile_id[1]*ts[1] + inner_raw[1] - iexch*ng[1];
                int k  = tile_id[2]*ts[2] + inner_raw[2] - iexch*ng[2];
                int lb = outer_raw[2];
                bool valid = true;
                auto nxx = nx;
                bool d = is_3d;
                const auto ngg  = ng;
                if constexpr (!is_3d)
                {
                    valid = valid && i <  nxx[0] + iexch*ngg[0];
                    valid = valid && j <  nxx[1] + iexch*ngg[1];
                    valid = valid && k <  nxx[2] + iexch*ngg[2];
                    valid = valid && i >= 0 - iexch*ngg[0];
                    valid = valid && j >= 0 - iexch*ngg[1];
                    valid = valid && k >= 0 - iexch*ngg[2];
                }
                index_type index(i, j, k, lb);
                if (valid)
                {
                    const val_t data = invoke_at(g_image, d_image, index, kernel);
                    d_image.set_elem(index, data);
                }
            });
        };
        dispatch::execute(outer_range, loop_load, kpool);
    }
    

    template <class array_t, class kernel_t>
    void fill_array(array_t& arr, const kernel_t& kernel, const grid::exchange_inclusion_e& exchange_policy = grid::include_exchanges)
    {
        transform_inplace(arr, kernel, exchange_policy);
    }
}