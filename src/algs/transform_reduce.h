#pragma once

#include "grid/grid.h"
#include "dispatch/support_of.h"
#include "dispatch/device_type.h"
#include "dispatch/execute.h"
#include "omni/omni.h"
#include "core/utils.h"
#include "algs/destructive_reduce.h"

namespace spade::algs
{
    const static struct max_t
    {
        template <typename data_t>
        _sp_hybrid auto operator() (const data_t& x0, const data_t& x1) const
        {
            return utils::max(x0, x1);
        }
    } max;
    
    const static struct sum_t
    {
        template <typename data_t>
        _sp_hybrid auto operator() (const data_t& x0, const data_t& x1) const
        {
            return x0 + x1;
        }
    } sum;
    
    template <typename func_t, typename operation_t, typename data_t, typename device_t>
    struct reduction_t
    {
        device_t idevice;
        const func_t func;
        const operation_t binary_op;
        mutable device::auto_vector<data_t, device_t> buffer;
        device_t device() const { return idevice; }
        void set_buf_size(const std::size_t n_elems) const { buffer.resize(n_elems); }
    };
    
    template <typename array_t, typename func_t, typename oper_t>
    inline auto make_reduction(const array_t& array, const func_t& func, const oper_t& oper)
    {
        const auto o_kern = omni::to_omni<array_t::centering_type()>(func, array);
        using data_t      = typename utils::remove_all<decltype(o_kern)>::type::output_type;
        using device_t    = typename array_t::device_type;
        device::auto_vector<data_t, device_t> tmp_buffer;
        return reduction_t<func_t, oper_t, data_t, device_t>{array.device(), func, oper, std::move(tmp_buffer)};
    }
    
    template <grid::multiblock_array array_t, typename func_t, typename operation_t, typename data_t>
    data_t transform_reduce(
        const array_t& array,
        const reduction_t<func_t, operation_t, data_t, typename array_t::device_type>& reduc)
    {
        const auto ctr       = array_t::centering_type();
        auto kernel          = omni::to_omni<ctr>(reduc.func, array);
        
        using kernel_t       = decltype(kernel);
        using output_t       = typename kernel_t::output_type;
        
        // This condition should be guaranteed if make_reduction is used
        static_assert(std::same_as<output_t, data_t>, "return type of reduction kernel must match buffer type!");        
        
        const auto& grid = array.get_grid();
        
        int pow_x = 2;
        int pow_y = 2;
        int pow_z = 2;
        
        int total_pow = pow_x + pow_y + pow_z;
        
        int nthr_x = 1 << pow_x;
        int nthr_y = 1 << pow_y;
        int nthr_z = 1 << pow_z;
        
        constexpr int dim = grid.dim();
        if (dim < 3) nthr_z = 1;
        
        int nblkx = utils::i_div_up(grid.get_num_cells(0), nthr_x);
        int nblky = utils::i_div_up(grid.get_num_cells(1), nthr_y);
        int nblkz = utils::i_div_up(grid.get_num_cells(2), nthr_z);
        
        int nblks = nblkx*nblky*nblkz;
        int nthrd = nthr_x*nthr_y*nthr_z;
        
        //Compute number of elements of the reduction buffer
        const std::size_t buffer_size = grid.get_num_local_blocks()*nblks;
        reduc.set_buf_size(buffer_size);
        
        auto buf_img = utils::make_vec_image(reduc.buffer);
        auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<output_t>(nthrd), dispatch::shmem::vec<bool>(nthrd));
        spade::dispatch::kernel_threads_t kpool(dispatch::ranges::make_range(0, nthr_x, 0, nthr_y, 0, nthr_z), array.device());
        
        int nx = grid.get_num_cells(0);
        int ny = grid.get_num_cells(1);
        int nz = grid.get_num_cells(2);
        
        int nx_max = utils::max(nx, ny, nz);
        
        if ((nx%2 != 0) || (nz%2 != 0) || ((nz%2 != 0) && (dim == 3)))
        {
            throw except::sp_exception("Attempted reduce operation on an array with odd number of cells per block, which is unsupported");
        }
        
        auto range = dispatch::ranges::make_range(0, nblkx, 0, nblky, 0, nblkz*grid.get_num_local_blocks());
        
        using index_t      = decltype(range)::index_type;
        using shmem_t      = decltype(k_shmem);
        using threads_t    = decltype(kpool);
        using inner_idx_t  = typename threads_t::index_type;
        
        const auto arr_img  = array.image();
        const auto grid_img = grid.image(array.device());
        
        int nbl = grid.get_num_local_blocks();
        
        auto loop = [=] _sp_hybrid (const index_t& idx, const threads_t& threads, shmem_t& shmem) mutable
        {
            //First start by filling the shmem vector with the results of the kernel invocation
            int i_blk  = idx[0];
            int j_blk  = idx[1];
            int k_blk  = idx[2] % nblkz;
            int lb_loc = (idx[2] - k_blk)/nblkz;
            
            auto& block_result_vec = shmem[0_c];
            auto& block_mask_vec   = shmem[1_c];
            
            const auto is_valid = [=](const grid::cell_idx_t& i_cell)
            {
                return (i_cell.i() < nx) && (i_cell.j() < ny) && (i_cell.k() < nz);
            };
            
            const auto block_shmem_idx = [=](const auto& ii)
            {
                return ii[0] + ii[1]*nthr_x + ii[2]*nthr_x*nthr_y;
            };
            
            threads.exec([&](const inner_idx_t& blk_idx)
            {
                // grid::cell_idx_t i_cell(i_blk + blk_idx[0]*nthr_x, j_blk + blk_idx[1]*nthr_y, k_blk + blk_idx[2]*nthr_z, lb_loc);
                grid::cell_idx_t i_cell(blk_idx[0] + i_blk*nthr_x, blk_idx[1] + j_blk*nthr_y, blk_idx[2] + k_blk*nthr_z, lb_loc);
                int shmem_idx = block_shmem_idx(blk_idx);
                if (is_valid(i_cell))
                {
                    auto val = invoke_at(grid_img, arr_img, i_cell, kernel);
                    block_mask_vec[shmem_idx] = true;
                    block_result_vec[shmem_idx] = val;
                }
                else
                {
                    block_mask_vec[shmem_idx] = false;
                }
            });
            
            threads.sync();
            //By this point, all of the results of the kernel are stored inside block_result_vec, so we do a 1-D reduction with the mask in mind
            for (int pitch_pow = 0; pitch_pow < total_pow; ++pitch_pow)
            {
                int pitch = 1 << pitch_pow;
                int mod   = pitch << 1;
                threads.exec([&](const inner_idx_t& blk_idx)
                {
                    int shmem_idx = block_shmem_idx(blk_idx);
                    if (shmem_idx % mod == 0)
                    {
                        int shmem_idx_2 = shmem_idx + pitch;
                        if (block_mask_vec[shmem_idx_2] && block_mask_vec[shmem_idx])
                        {
                            block_result_vec[shmem_idx] = reduc.binary_op(block_result_vec[shmem_idx], block_result_vec[shmem_idx_2]);
                        }
                    }
                });
            }
            
            if (threads.isroot())
            {
                std::size_t vec_idx = i_blk + nblkx*(j_blk + nblky*(k_blk + nblkz*lb_loc));
                buf_img[vec_idx] = block_result_vec[0];
            }
        };
        
        dispatch::execute(range, loop, kpool, k_shmem);
        
        //Compute result on this rank
        const auto local_result = destructive_reduce(reduc.buffer, reduc.binary_op);
        const auto& group = grid.group();
        return group.reduce(local_result, reduc.binary_op);
    }
}