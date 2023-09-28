#pragma once

#include "core/coord_system.h"

namespace spade::grid
{
    template <
        typename uelem_t,
        coords::coordinate_system coord_t,
        template <typename> typename container_tt
        >
    struct unstructured_geometry_t
    {
        //todo
        using float_t = coord_t::coord_type;
        using point_t = coords::point_t<float_t>;
        using uint_t  = std::size_t;
        container_tt<point_t> points;
    };
}