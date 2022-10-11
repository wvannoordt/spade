#pragma once
#include "spade.h"

using grid_t = spade::grid::cartesian_grid_t
<
    spade::coords::identity<double>,
    spade::parallel::mpi_t,
    3
>;

using array_t = spade::grid::cell_array
<
    grid_t,
    double,
    spade::dims::singleton_dim,
    std::vector<typename spade::grid::detail::get_fundamental_type<double>::type>
>;

void diffkernel(int i, int j, int k, int lb, array_t& rhs, const array_t& q);