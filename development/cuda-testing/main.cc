#include <chrono>
#include <print.h>
//#include <spade.h>

using real_t = double;

__host__ __device__ real_t twice(const real_t r) {return 2.0*r;} 



template <typename thing_t> __global__ void kern(thing_t val)
{
    auto pp = twice(val);
}

int main(int argc, char** argv)
{
    kern<<<2,2>>>(0.1);
    print(cudaGetErrorString(cudaGetLastError()));
    return 0;
}
