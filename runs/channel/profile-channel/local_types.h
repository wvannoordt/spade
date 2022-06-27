#pragma once

#include "spade.h"
#include "prof_t.h"

typedef double real_t;
typedef spade::ctrs::array<real_t, 3> v3d;
typedef spade::ctrs::array<real_t, 3> v3r;
typedef spade::ctrs::array<real_t, 6> v6r;
typedef spade::ctrs::array<int,    3> v3i;
typedef spade::ctrs::array<int,    4> v4i;
typedef spade::ctrs::array<spade::grid::cell_t<int>, 4> v4c;
typedef spade::ctrs::array<spade::grid::face_t<int>, 5> v5f;
typedef spade::fluid_state::prim_t<real_t> prim_t;
typedef spade::fluid_state::cons_t<real_t> cons_t;
typedef postprocessing::prof_t<real_t> profr_t;
typedef spade::linear_algebra::dense_mat<real_t, 3> m3r;