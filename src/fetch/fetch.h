#pragma once

#include "fetch/data_types.h"
#include "fetch/collections.h"
#include "fetch/stencils.h"
#include "fetch/fetch_types.h"
#include "fetch/retrieval.h"

namespace spade::fetch
{
    template <grid::multiblock_grid grid_t, grid::multiblock_array array_t, typename flux_in_t>
    void get_face_data(const grid_t& grid, const array_t& prims, const grid::face_idx_t& iface, flux_in_t& flux_input)
    {
        detail::get_face_stencil_data(grid, prims, iface, flux_input.cell_data);
        detail::get_single_face_data(grid, prims, iface, flux_input.face_data);
    }
    
    template <grid::multiblock_grid grid_t, grid::multiblock_array array_t, typename cell_in_t>
    void get_cell_data(const grid_t& grid, const array_t& prims, const grid::cell_idx_t& icell, cell_in_t& input)
    {
        detail::get_cell_stencil_data(grid, prims, icell, input.cell_data);
    }
}