#pragma once

#include <concepts>

#include "grid/grid.h"
#include "omni/omni.h"

namespace spade::time_integration
{
    template <typename rhs_arr_t, typename cons_arr_t, typename limiter_t>
    static void advance_limiter(rhs_arr_t& dsol, const cons_arr_t& sol, const limiter_t& limiter_func)
    {
		// Setup image of arrays so the GPU can access it
        const auto cons_arr = sol.image();
        auto dcons_arr      = dsol.image();

		// Get looping range
        auto grid_range = dispatch::support_of(sol, grid::exclude_exchanges);
        auto loop = [=] _sp_hybrid (const grid::cell_idx_t& icell) mutable
        {
			// Get data from image
			const auto cons_val  = cons_arr.get_elem(icell);
			auto dcons_val = dcons_arr.get_elem(icell);
			
			// Run limiter on this grid point
			auto update = limiter_func(dcons_val, cons_val);

			// Set new update in array
            dcons_arr.set_elem(icell, update);
        };
		// Call the GPU to evaluate this
        dispatch::execute(grid_range, loop);
    }
}
