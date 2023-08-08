#pragma once

#ifdef __NVCC__
#define _sp_cuda 1
#else
#define _sp_cuda 0
#endif

#if (_sp_cuda)

#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>

#define _sp_device __device__
#define _sp_hybrid __host__ __device__

//This has been defined as _sp_hybrid instead of _sp_device because static analysis of the closure type is not
//correct for just __device__ lambdas.
//See https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#notes-on-host-device-lambdas
//Section 14.7.3 and 14.7.2.14
#define _sp_lambda [=] _sp_hybrid

#define _sp_cu_guard(mycode) mycode
#else

#define _sp_device
#define _sp_hybrid
#define _sp_lambda [=]
#define _sp_cu_guard(mycode) 

#endif