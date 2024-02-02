#pragma once

#include "core/vec_image.h"
#include "dispatch/execute.h"

namespace spade::algs
{
    //Note that this should actually be the return type of the unary operation!
    template <typename container_t, typename binary_t, typename unary_t>
    inline typename container_t::value_type
    destructive_reduce(container_t& vec, const binary_t& binary_op, const unary_t& unary_op)
    {
        using data_t = typename container_t::value_type;
        int threads_pow = 3;
        std::size_t num_threads = (1 << threads_pow);
        
        auto vec_img = utils::make_vec_image(vec);
        auto device  = device::from_vector(vec);
        
        std::size_t extent_size = vec_img.size();
        while (extent_size > 1)
        {
            auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<data_t>(num_threads));
            auto range   = dispatch::ranges::make_range(0, utils::i_div_up(extent_size, num_threads));
            spade::dispatch::kernel_threads_t kpool(dispatch::ranges::make_range(0, num_threads), device);
            using shmem_t   = decltype(k_shmem);
            using threads_t = decltype(kpool);
            auto loop = [=] _sp_hybrid (const std::size_t& major, const threads_t& threads, shmem_t& shmem) mutable
            {
                auto& shvec = shmem[0_c];
                
                threads.exec([&](const std::size_t& minor)
                {
                    std::size_t global_idx = minor + major*num_threads;
                    if (global_idx < extent_size)
                    {
                        shvec[minor] = unary_op(vec_img[global_idx]);
                    }
                });
                
                for (int pitch_pow = 0; pitch_pow < threads_pow; ++pitch_pow)
                {
                    int pitch = 1 << pitch_pow;
                    int mod   = pitch << 1;
                    threads.exec([&](const std::size_t& minor)
                    {
                        std::size_t global_idx  = minor + major*num_threads;
                        std::size_t global_idx2 = global_idx + pitch;
                        if (((global_idx % mod) == 0) && (global_idx < extent_size) && ((global_idx2 < extent_size)))
                        {
                            shvec[minor] = binary_op(shvec[minor], shvec[minor + pitch]);
                        }
                    });
                    
                    threads.sync();
                }
                
                if (threads.isroot())
                {
                    vec_img[major] = shvec[0];
                }
            };
            
            dispatch::execute(range, loop, kpool, k_shmem);
            extent_size = utils::i_div_up(extent_size, num_threads);
        }
        return device::inspect_at(vec, 0);
    }
    
    template <typename container_t, typename binary_t>
    inline typename container_t::value_type
    destructive_reduce(container_t& vec, const binary_t& binary_op)
    {
        return destructive_reduce(vec, binary_op, [] _sp_hybrid (const typename container_t::value_type& in) { return in; });
    }
}