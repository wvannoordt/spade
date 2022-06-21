#pragma once

#include "cvdf.h"
#include "prof_t.h"

typedef double real_t;
typedef cvdf::ctrs::array<real_t, 3> v3d;
typedef cvdf::ctrs::array<real_t, 3> v3r;
typedef cvdf::ctrs::array<real_t, 6> v6r;
typedef cvdf::ctrs::array<int,    3> v3i;
typedef cvdf::ctrs::array<int,    4> v4i;
typedef cvdf::ctrs::array<cvdf::grid::cell_t<int>, 4> v4c;
typedef cvdf::ctrs::array<cvdf::grid::face_t<int>, 5> v5f;
typedef cvdf::fluid_state::prim_t<real_t> prim_t;
typedef cvdf::fluid_state::cons_t<real_t> cons_t;
typedef postprocessing::prof_t<real_t> profr_t;
typedef cvdf::linear_algebra::dense_mat<real_t, 3> m3r;