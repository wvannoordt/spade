#pragma once

#include <cmath>

#include "core/config.h"
#include "core/typedef.h"
#include "core/consts.h"
#include "core/print.h"
#include "core/range.h"
#include "core/typedef.h"
#include "core/grid.h"
#include "core/io.h"
#include "core/parallel.h"
#include "core/algs.h"
#include "core/timing.h"
#include "core/reduce_ops.h"
#include "core/linear_algebra.h"
#include "core/time_integration.h"
#include "core/block_config.h"

#include "navier-stokes/fluid_state.h"
#include "navier-stokes/viscous_laws.h"
#include "navier-stokes/convective.h"
#include "navier-stokes/viscous.h"
#include "navier-stokes/pde_algs.h"
#include "navier-stokes/navier_stokes_mms.h"