#include <chrono>
#include <spade.h>

using real_t = double;

_spade_hybrid_ real_t twice(const real_t r) {return 2.0*r;} 

#if (_spade_cuda_)

template <typename thing_t> __global__ void kern(thing_t val)
{
    auto pp = twice(val);
}

#endif

int main(int argc, char** argv)
{
    kern<<<2,2>>>(0.1);
    print(cudaGetErrorString(cudaGetLastError()));
    return 0;
}
