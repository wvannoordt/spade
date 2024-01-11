#include "spade.h"
#include "dispatch/execute2.h"
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
    auto bnd_outer = spade::utils::make_bounds(0, 500);
    auto bnd_inner = spade::utils::make_bounds(0, 32);
    auto range     = spade::dispatch::ranges::make_range(bnd_outer);
    spade::dispatch::kernel_threads_t kpool(range, device);
    
    using threads_t = decltype(kpool);
    const auto loop = [=] _sp_hybrid (const std::size_t& idx, const threads_t& threads)
    {
        threads.exec([&](const auto& i){ printf("i: %d", i); });
    };
    
    
    
    spade::dispatch::proto::execute(range, loop, kpool);
    
    
    return 0;
}
