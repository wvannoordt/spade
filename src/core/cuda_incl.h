#pragma once

#ifdef __NVCC__
#define _spade_cuda_ 1
#else
#define _spade_cuda_ 0
#endif

#if (_spade_cuda_)

#define _spade_device_ __device__
#define _spade_hybrid_ __host__ __device__

#else

#define _spade_device_
#define _spade_hybrid_

#endif