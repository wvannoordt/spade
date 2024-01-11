#pragma once

#include "core/cuda_incl.h"
#include "dispatch/ranges/kernel_config.h"
#include "dispatch/ranges/grid_idx_range.h"
#include "dispatch/ranges/linear_range.h"

namespace spade::dispatch
{
    
    template <typename idx_t, typename range_t> struct kernel_params_t{};
    
    template <typename idx_t, typename dev_t>
    struct kernel_params_t<idx_t, ranges::grid_idx_range_t<idx_t, dev_t>>
    {
        using range_type = ranges::grid_idx_range_t<idx_t, dev_t>;
        using range_t      = range_type;
        using index_type   = idx_t;
        using triplet_type = ranges::d3_t;
        
        constexpr static int blocksize = 4;
        triplet_type grid_size;
        triplet_type block_size;
        range_type irange;
        
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
    
    template <typename idx_t, typename dev_t>
    struct kernel_params_t<idx_t, ranges::linear_range_t<idx_t, dev_t>>
    {
        using range_type   = ranges::linear_range_t<idx_t, dev_t>;
        using range_t      = range_type;
        using index_type   = typename range_t::index_type;
        using triplet_type = ranges::d3_t;
        
        constexpr static int blocksize = 32;
        triplet_type grid_size;
        triplet_type block_size;
        range_type irange;
        
        _sp_hybrid index_type get_index(const triplet_type& thr_idx, const triplet_type& blo_idx, const triplet_type& blo_dim, const triplet_type& grd_dim) const
        {
            index_type output;
            output[0] = irange.lower + thr_idx.x + blo_idx.x*blo_dim.x;
            return output;
        }
        
        _sp_hybrid bool is_valid(const index_type& i) const
        {
            return (i[0] < irange.upper) && (i[0] >= irange.lower);
        }
    };
    
    template <typename idx_t, typename dev_t>
    struct kernel_params_t<idx_t, ranges::bilinear_range_t<idx_t, dev_t>>
    {
        using range_type   = ranges::bilinear_range_t<idx_t, dev_t>;
        using range_t      = range_type;
        using index_type   = typename range_t::index_type;
        using triplet_type = ranges::d3_t;
        
        constexpr static int blocksize = 16;
        triplet_type grid_size;
        triplet_type block_size;
        range_type irange;
        
        _sp_hybrid index_type get_index(const triplet_type& thr_idx, const triplet_type& blo_idx, const triplet_type& blo_dim, const triplet_type& grd_dim) const
        {
            index_type output;
            output[0] = irange.i0 + thr_idx.x + blo_idx.x*blo_dim.x;
            output[1] = irange.j0 + thr_idx.y + blo_idx.y*blo_dim.y;
            return output;
        }
        
        _sp_hybrid bool is_valid(const index_type& i) const
        {
            return (i[0] < irange.i1) && (i[0] >= irange.i0) && (i[1] < irange.j1) && (i[1] >= irange.j0);
        }
    };
    
    template <grid::grid_index index_t, typename r_device_t>
    static auto get_launch_params(const ranges::grid_idx_range_t<index_t, r_device_t>& range_in)
    {
        using range_t = ranges::grid_idx_range_t<index_t, r_device_t>;
        using index_type = typename range_t::index_type;
        using triplet_type = typename kernel_params_t<index_type, range_t>::triplet_type;
        
        
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
    
    template <typename index_t, typename r_device_t>
    static auto get_launch_params(const ranges::linear_range_t<index_t, r_device_t>& range_in)
    {
        using range_t = ranges::linear_range_t<index_t, r_device_t>;
        using index_type = typename range_t::index_type;
        using triplet_type = typename kernel_params_t<index_t, range_t>::triplet_type;
        
        constexpr int blocksize = kernel_params_t<index_t, range_t>::blocksize;
        
        triplet_type gsize;
        gsize.x =  1 + ((range_in.upper - range_in.lower - 1) / blocksize);
        gsize.y =  1;
        gsize.z =  1;
        
        triplet_type bsize;
        bsize.x = blocksize;
        bsize.y = 1;
        bsize.z = 1;        
        
        return kernel_params_t<index_t, range_t>{gsize, bsize, range_in};
    }
    
    template <typename index_t, typename r_device_t>
    static auto get_launch_params(const ranges::bilinear_range_t<index_t, r_device_t>& range_in)
    {
        using range_t = ranges::bilinear_range_t<index_t, r_device_t>;
        using index_type = typename range_t::index_type;
        using triplet_type = typename kernel_params_t<index_t, range_t>::triplet_type;
        
        constexpr int blocksize = kernel_params_t<index_t, range_t>::blocksize;
        
        triplet_type gsize;
        gsize.x =  1 + ((range_in.i1 - range_in.i0 - 1) / blocksize);
        gsize.y =  1 + ((range_in.j1 - range_in.j0 - 1) / blocksize);
        gsize.z =  1;
        
        triplet_type bsize;
        bsize.x = blocksize;
        bsize.y = blocksize;
        bsize.z = 1;        
        
        return kernel_params_t<index_t, range_t>{gsize, bsize, range_in};
    }
}