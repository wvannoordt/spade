#pragma once

#include "dispatch/shared_vector.h"
#include "grid/grid_index_types.h"

namespace spade::ibm
{
    template <typename float_t>
    struct ghost_list_t
    {
        using idx_t = grid::cell_idx_t;
        using pnt_t = coords::point_t<float_t>;
        using vec_t = ctrs::array<float_t, pnt_t::size()>;
        
        //Stores the points locally
        device::shared_vector<idx_t>       indices;
        device::shared_vector<pnt_t>       boundary_points;
        device::shared_vector<vec_t>       boundary_normals;
        device::shared_vector<pnt_t>       closest_points;
        device::shared_vector<vec_t>       closest_normals;
        device::shared_vector<int>         directions;
        device::shared_vector<int>         signs;
    };
}