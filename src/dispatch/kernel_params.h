#pragma once

#include "core/cuda_incl.h"

namespace spade::dispatch
{
    
    struct dim3_wrapper_t
    {
        int x, y, z;
    };
    
#if (_sp_cuda)
    using d3_t = dim3;
#else
    using d3_t = dim3_wrapper_t;
#endif
    
    template <typename idx_t, typename range_t>
    struct kernel_params_t
    {
        using index_type   = idx_t;
        using triplet_type = d3_t;
        
        constexpr static int blocksize = 8;
        triplet_type grid_size;
        triplet_type block_size;
        range_t irange;
        
        _sp_hybrid index_type get_index(const triplet_type& thr_idx, const triplet_type& blo_idx, const triplet_type& blo_dim, const triplet_type& grd_dim) const
        {
            index_type output;
            auto raw_x   = thr_idx.x + blo_idx.x*blo_dim.x;
            const int nx = irange.upper.i()-irange.lower.i();
            const int sx = blocksize*(1+(nx-1)/blocksize);
            auto raw_lb = raw_x/sx;
            output.i()  = irange.lower.i()  + raw_x - raw_lb*sx;
            output.j()  = irange.lower.j()  + thr_idx.y + blo_idx.y*blo_dim.y;
            output.k()  = irange.lower.k()  + thr_idx.z + blo_idx.z*blo_dim.z;
            output.lb() = irange.lower.lb() + raw_lb;
            return output;
        }
        
        _sp_hybrid bool is_valid(const index_type& i) const
        {
            bool output = true;
            algs::static_for<0, index_type::size()>([&](const auto& d)
            {
                output = output && (i[d.value] >= irange.lower[d.value] && i[d.value] < irange.upper[d.value]);
            });
            return output;
        }
    };
    
    template <typename range_t>
    static auto get_launch_params(const range_t& range_in)
    {
        using index_type = typename range_t::index_type;
        using triplet_type = kernel_params_t<index_type, range_t>::triplet_type;
        
        
        int ni = range_in.upper.i() - range_in.lower.i();
        int nj = range_in.upper.j() - range_in.lower.j();
        int nk = range_in.upper.k() - range_in.lower.k();
        
        constexpr int blocksize = kernel_params_t<index_type, range_t>::blocksize;
        
        triplet_type gsize;
        gsize.x =  1+((ni-1)/blocksize);
        gsize.x *= (range_in.upper.lb() - range_in.lower.lb());
        gsize.y =  1+((nj-1)/blocksize);
        gsize.z =  1+((nk-1)/blocksize);
        
        triplet_type bsize;
        bsize.x = blocksize;
        bsize.y = blocksize;
        bsize.z = blocksize;
        if (range_in.upper.k() - range_in.lower.k() == 1) bsize.z = 1;
        
        return kernel_params_t<index_type, range_t>{gsize, bsize, range_in};
    }
}