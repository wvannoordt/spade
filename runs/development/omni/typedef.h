#pragma once

using real_t = double;
using v3d    = spade::ctrs::array<real_t, 3>;
using v3i    = spade::ctrs::array<int,    3>;
using v4i    = spade::ctrs::array<int,    4>;
using flux_t = spade::fluid_state::flux_t<real_t>;
using prim_t = spade::fluid_state::prim_t<real_t>;
using cons_t = spade::fluid_state::cons_t<real_t>;
