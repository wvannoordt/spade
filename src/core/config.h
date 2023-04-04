#pragma once

#if(__cplusplus>=202002)
#define C20_G11 1
#ifndef USE_SOURCE_LOCATION
#define USE_SOURCE_LOCATION 0
#endif
#else
#define C20_G11 0
#define USE_SOURCE_LOCATION 0
#endif

#include "core/c20.h"