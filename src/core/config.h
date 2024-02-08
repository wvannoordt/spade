#pragma once

#if(__cplusplus>=202002)
#define C20_G11 1
#else
#define C20_G11 0
#endif

#include "cuda_incl.h"

#if (_sp_cuda)
#   if (CUDART_VERSION >= 12020)
#       define USE_SOURCE_LOCATION 1
#   else
#       define USE_SOURCE_LOCATION 0
#   endif
#else
#   define USE_SOURCE_LOCATION C20_G11
#endif

#include "core/c20.h"