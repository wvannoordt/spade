#include "spade.h"
using real_t = double;

int main(int argc, char** argv)
{
    const auto device = spade::device::gpu;
    
    // //Basic usage
    // auto bnd    = spade::utils::make_bounds(0, 5, 0, 5);
    // auto range  = spade::dispatch::ranges::make_range(bnd);
    
    // using idx_t = decltype(range)::index_type;
    // const auto loop = [=] _sp_hybrid (const idx_t& idx)
    // {
    //     //do something wiht idx.
    // };
    
    // spade::dispatch::execute(range, loop);
    
    
    //"composite" usage, to enforce cooperation between block threads
    //Note that you need to have BOTH the threads and the shared data arguments
    auto bnd_outer = spade::utils::make_bounds(0, 6);
    auto bnd_inner = spade::utils::make_bounds(0, 5);
    auto o_range   = spade::dispatch::ranges::make_range(bnd_outer);
    auto i_range   = spade::dispatch::ranges::make_range(bnd_inner);
    spade::dispatch::kernel_threads_t kpool(i_range, device);
    //   0  1  2  3  4
    //0: [x  x  x  x  x]
    //1: x  x  x  x  x
    //2: x  x  x  x  x
    // x  x  x  x  x
    // x  x  x  x  x
    // x  x  x  x  x
    
    
    using threads_t = decltype(kpool);
    const auto loop = [=] _sp_hybrid (const std::size_t& idx, const threads_t& threads)
    {
        threads.exec([&](const auto& i){ printf("i, idx: %d, %d\n", i, int(idx)); });
        threads.sync();
        auto result0 = threads.sum([&](const auto& i) { return real_t(i + 100*idx); });
        auto result1 = threads.max([&](const auto& i) { return real_t(i + 100*idx); });
        auto result2 = threads.min([&](const auto& i) { return real_t(i + 100*idx); });
        threads.exec([&](const auto& i){ printf("%f, %f, %f\n", result0, result1, result2); });
    };
    
    spade::dispatch::execute(o_range, loop, kpool);
    
    auto shmem = spade::dispatch::shmem::make_shmem(
        spade::dispatch::shmem::vec<real_t>(bnd_inner.volume()),
        spade::dispatch::shmem::vec<int>   (bnd_inner.volume()),
        spade::dispatch::shmem::vec<char>  (bnd_inner.volume()));
    
    using shmem_t = decltype(shmem);
    
    const auto loop2 = [=] _sp_hybrid (const std::size_t& idx, const threads_t& threads, shmem_t& shmem)
    {
        auto& v_r = shmem[0_c];
        auto& v_i = shmem[1_c];
        auto& v_c = shmem[2_c];
        
        threads.exec([&](const auto i) { v_r[i] = real_t(i)+real_t(0.25); });
        threads.exec([&](const auto i) { v_i[i] = -i; });
        threads.exec([&](const auto i) { v_c[i] = i+1; });
        threads.sync();
        threads.exec([&](const auto i) { printf("%f\n", v_r[i]);});
        threads.sync();
        threads.exec([&](const auto i) { printf("%d\n", v_i[i]);});
        threads.sync();
        threads.exec([&](const auto i) { printf("%d\n", v_c[i]);});
        threads.sync();
        threads.exec([&](const auto i) { printf("%d\n", v_i[i]);});
        threads.sync();
        threads.exec([&](const auto i) { printf("%f\n", v_r[i]);});
        threads.sync();
    };
    
    spade::dispatch::execute(o_range, loop2, kpool, shmem);
    
    return 0;
}
