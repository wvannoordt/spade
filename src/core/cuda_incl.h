#pragma once

#ifdef __NVCC__
#define _spade_cuda_ 1
#else
#define _spade_cuda_ 0
#endif

#if (_spade_cuda_)

#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>

#define _sp_device __device__
#define _sp_hybrid __host__ __device__
#define _sp_lambda [=] _sp_hybrid
#else

#define _sp_device
#define _sp_hybrid
#define _sp_lambda [=]

#endif