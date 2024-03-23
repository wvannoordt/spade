#pragma once

#include "grid/grid.h"
#include "dispatch/support_of.h"
#include "dispatch/device_type.h"
#include "dispatch/execute.h"
#include "omni/omni.h"
#include "core/invoke.h"

namespace spade::algs
{
    /*
    template <grid::multiblock_array array_t, class callable_t>
    void transform_inplace(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy=grid::exclude_exchanges)
    {
        const auto ctr       = array_t::centering_type();
        auto kernel          = omni::to_omni<ctr>(func, arr);
        auto var_range       = dispatch::support_of(arr, exchange_policy);
        using index_type     = typename decltype(var_range)::index_type;
        auto d_image         = arr.image();
        const auto g_image   = arr.get_grid().image(arr.device());
        
        using val_t = typename array_t::alias_type;
        auto loop_load = [=] _sp_hybrid (const index_type& index) mutable
        {
            const val_t data = invoke_at(g_image, d_image, index, kernel);
            d_image.set_elem(index, data);
        };
        dispatch::execute(var_range, loop_load);
    }
    */
    
    template <grid::multiblock_array array_t, class callable_t>
    void transform_inplace(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy=grid::exclude_exchanges)
    {
        const auto ctr       = array_t::centering_type();
        auto kernel          = omni::to_omni<ctr>(func, arr);
        auto var_range       = dispatch::support_of(arr, exchange_policy);
        using index_type     = typename decltype(var_range)::index_type;
        auto d_image         = arr.image();
        const auto g_image   = arr.get_grid().image(arr.device());
        
        using val_t = typename array_t::alias_type;
        auto loop_load = [=] _sp_hybrid (const index_type& index) mutable
        {
            const val_t data = invoke_at(g_image, d_image, index, kernel);
            d_image.set_elem(index, data);
        };
        dispatch::execute(var_range, loop_load);
    }

    template <class array_t, class kernel_t>
    void fill_array(array_t& arr, const kernel_t& kernel, const grid::exchange_inclusion_e& exchange_policy = grid::include_exchanges)
    {
        transform_inplace(arr, kernel, exchange_policy);
    }
}