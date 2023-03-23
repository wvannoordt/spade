#pragma once

#include <type_traits>

#include "core/config.h"
#include "core/attribs.h"
#include "core/grid.h"
#include "core/reduce_ops.h"
#include "core/md_loop.h"
#include "core/invoke.h"

namespace spade::algs
{
    template <typename array_t, typename kernel_t>
    void block_loop(array_t& array, const int& lb, const kernel_t& kernel, const grid::exchange_inclusion_e& exchange_policy)
    {
        auto bound_box = array.get_grid_index_bounding_box(exchange_policy);
        const int lb_dim = array_t::index_type::lb_idx;
        
        //modify the bounding box to only loop over this one block
        bound_box.min(lb_dim) = (int)(lb);
        bound_box.max(lb_dim) = (int)(lb+1);
        md_loop<int, 4, kernel_t, typename array_t::index_type>(bound_box, kernel);
    }

    template <grid::multiblock_array array_t, class callable_t, reduce_ops::reduce_operation<typename array_t::value_type> reduce_t>
    auto transform_reduce(
        const array_t& arr,
        const callable_t& func,
        reduce_t& reduce_oper,
        const grid::exchange_inclusion_e& exchange_policy=grid::exclude_exchanges)
    {
        const auto& grid = arr.get_grid();
        const auto nlb = grid.get_num_local_blocks();
        const typename array_t::index_type ii = 0;
        auto init_data = arr.get_elem(ii);
        reduce_oper.init(func(init_data));
        for (auto lb: range(0, nlb))
        {
            const auto loop_func = [&](const auto& elem) -> void
            {
                const auto data = arr.get_elem(elem);
                reduce_oper.reduce_elem(func(data));
            };
            block_loop(arr, lb, loop_func, exchange_policy);
        }
        return grid.group().reduce(reduce_oper.value,reduce_oper.equiv_par_op());
    }
    
    template <grid::multiblock_array array_t, class callable_t>
    auto& transform_inplace(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy=grid::exclude_exchanges)
    {
        const auto& grid = arr.get_grid();
        const auto nlb = grid.get_num_local_blocks();

        const auto kernel = omni::to_omni(func);
        
        // consider fundamentally separating the block dimension with the ijk dimensions!
        for (auto lb: range(0, nlb))
        {
            const auto loop_func = [&](const auto& elem) -> void
            {
                const auto data = invoke_at(func, arr, elem);
                arr.set_elem(elem, data);
            };
            block_loop(arr, lb, loop_func, exchange_policy);
        }
        return arr;
    }


    template <class array_t, class kernel_t>
    void fill_array(array_t& arr, const kernel_t& kernel, const grid::exchange_inclusion_e& exchange_policy=grid::include_exchanges)
    {
        transform_inplace(arr, kernel, exchange_policy);
    }
}
