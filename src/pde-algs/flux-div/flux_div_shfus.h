#pragma once

#include <concepts>

#include "grid/grid.h"
#include "omni/omni.h"

#include "pde-algs/pde_traits.h"
#include "pde-algs/flux-div/tags.h"

namespace spade::pde_algs
{
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t,
        typename traits_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    inline void flux_div_shfus(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func,        
        const traits_t& traits)
    {
        using real_type     = sol_arr_t::value_type;
        using alias_type    = sol_arr_t::alias_type;
        using omni_type     = typename flux_func_t::omni_type;
        using flux_type     = rhs_arr_t::alias_type;
        using grid_type     = typename sol_arr_t::grid_type;
        
        static_assert(std::same_as<typename grid_type::coord_sys_type, coords::identity<typename grid_type::coord_type>>, "flux divergence does not yet account for the jacobian!");
        
        using namespace sym::literals;
        
        // Decide at compile-time whether or not the flux divergence overwrites with zero before writing
        const auto& incr = algs::get_trait(traits, "pde_increment"_sym, increment);
        using incr_mode_t = typename utils::remove_all<decltype(incr)>::type;
        constexpr bool is_incr_mode = incr_mode_t::increment_mode;
        
        // Grid geometry is stored in the grid array
        const auto& grid    = prims.get_grid();
        
        // GPU-compatible "image" of the grid object
        const auto grid_img = grid.image(partition::local, prims.device());
        
        // GPU-compatible "image" of the flowfield and residual arrays
        const auto q_img    = prims.image();
        auto rhs_img        = rhs.image();
        
        // 2D or 3D simulation
        constexpr int dim = grid.dim();
        
        // Number of cells in a block
        auto nx     = grid.get_num_cells();
        
        // Number of exchange/halo cells in a block
        auto ngs    = prims.get_num_exchange();
        auto ntiles = nx;
        
        // Compute extents of the flux stencil
        constexpr int left_extent  = -static_math::moddiv<omni_type::template min_extent<0>, 2>::value;
        constexpr int right_extent =  static_math::moddiv<omni_type::template max_extent<0>, 2>::value + 1;
        
        // Compute number of required halo / exchange cells
        constexpr int ng           =  utils::max(left_extent, right_extent);
        constexpr int tile_size    =  4;
        constexpr int half_tile    =  tile_size / 2;
        
        // Compute number of tiles in each direction
        int total_tiles = 1;
        for (int d = 0; d < nx.size(); ++d)
        {
            ntiles[d]    = utils::i_div_up(nx[d], tile_size);
            total_tiles *= ntiles[d];
        }
        
        // Create data structure for the flux kernel
        using input_type = typename omni::stencil_data_t<omni_type, sol_arr_t>;
        
        // constexpr int thin_dir = 1;
        // constexpr int fuse_dir = 0;
        // constexpr int seq_dir  = 2;
        
        constexpr int thin_dir = 0;
        constexpr int fuse_dir = 2;
        constexpr int seq_dir  = 1;
        
        spade::ctrs::array<int, 3> irange_dims = tile_size;
        irange_dims[fuse_dir] *= 2;
        ntiles[fuse_dir] /= 2;
        ntiles[seq_dir]  /= 2;
        
        // Create a range for the individual tile size
        const auto tile_range = dispatch::ranges::make_range(0, irange_dims[0], 0, irange_dims[1], 0, irange_dims[2]);
        
        // Create collaborative threads object over the tile range
        dispatch::kernel_threads_t kpool(tile_range, prims.device());
        using threads_type = decltype(kpool);
        
        ctrs::array<int, 3> shsize;
        shsize[thin_dir] = tile_size;
        shsize[fuse_dir] = tile_size*2;
        shsize[seq_dir]  = tile_size*2;
        const std::size_t total_sh_vals = shsize[0]*shsize[1]*shsize[2];
        
        // Create shared memory array
        auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<flux_type>(total_sh_vals));
        using shmem_type = decltype(k_shmem);
        
        // Create a grid-range over the number of tiles in each block and the number of blocks
        const auto outer_range = dispatch::ranges::make_range(0, ntiles[0], 0, ntiles[1]*ntiles[2], 0, int(grid.get_num_local_blocks()));
        
        // Create a lambda for the bulk workload
        auto loop = [=] _sp_hybrid (const ctrs::array<int, 3>& outer_raw, const threads_type& threads, shmem_type& shmem) mutable
        {
            // Compute the tile index for this thread block
            int tile_id_1d = outer_raw[1];
            auto& shmem_vec = shmem[0_c];
            
            auto sh_rhs = utils::make_vec_image(shmem_vec, shsize);
            ctrs::array<int, 3> tile_id;
            tile_id[0]  = outer_raw[0];
            tile_id[1]  = tile_id_1d % ntiles[1];
            tile_id_1d -= tile_id[1];
            tile_id_1d /= ntiles[1];
            tile_id[2]  = tile_id_1d;

            input_type input0;
            input_type input;
            // Block index
            int lb = outer_raw[2];
            
            // 1 / (grid spacing) for differentiation
            const auto inv_dx_native = grid_img.get_inv_dx(lb);
            
            // Convert 1/deltaX into the necessary precision for this computation
            ctrs::array<real_type, dim> inv_dx;
            #pragma unroll
            for (int d = 0; d < dim; ++d) inv_dx[d] = inv_dx_native[d];
            
            threads.exec([&](const ctrs::array<int, 3>& inner_raw)
            {
                int dum0 = thin_dir;
                int dum1 = seq_dir;
                int dum2 = fuse_dir;
                #pragma unroll
                for (int t_id = 0; t_id < 2; ++t_id)
                {
                    auto ii = inner_raw;
                    ii[seq_dir] += tile_size*t_id;
                    flux_type& my_rhs = sh_rhs(ii[0], ii[1], ii[2]);
                    const auto nothing = rhs_img.size();
                    constexpr bool is_incr_modetmp = is_incr_mode;
                    if constexpr (is_incr_modetmp)
                    {
                        grid::cell_idx_t i_cell;
                        i_cell.lb()         = lb;
                        i_cell.i(thin_dir)  = tile_id[thin_dir ]*tile_size        + inner_raw[thin_dir];
                        i_cell.i(seq_dir )  = (2*tile_id[seq_dir]+t_id)*tile_size + inner_raw[seq_dir];
                        i_cell.i(fuse_dir)  = 2*tile_id[fuse_dir]*tile_size       + inner_raw[fuse_dir];
                        my_rhs = rhs_img.get_elem(i_cell);
                    }
                    else
                    {
                        my_rhs = real_type(0.0);
                    }
                    threads.sync();
                }
                
                const auto get_offsts = [&](const spade::ctrs::array<int, 3>& ii)
                {
                    spade::ctrs::array<int, 4> output;
                    bool is_seq_face  = ii[seq_dir] == (tile_size - 1);
                    bool is_fuse_face = ii[seq_dir] == (tile_size - 2);
                    bool is_big_face  = !(is_seq_face || is_fuse_face);
                    int flux_dir  = is_big_face ? thin_dir : (is_seq_face ? seq_dir : fuse_dir);
                    int big_dfuse = ii[fuse_dir];
                    int big_dseq  = ii[thin_dir] + tile_size * ii[seq_dir];
                    int lit_dthn  = ii[thin_dir];
                    int lit_dd    = ii[fuse_dir];
                    
                    int di_fuse_dir = is_big_face ? (big_dfuse)     : (is_seq_face  ? (lit_dd) : (2*tile_size - 1));
                    int di_seq_dir  = is_big_face ? (big_dseq)      : (is_fuse_face ? (lit_dd) : (2*tile_size - 1));
                    int di_thin_dir = is_big_face ? (tile_size - 1) : lit_dthn;
                    
                    output[fuse_dir] = di_fuse_dir;
                    output[seq_dir ] = di_seq_dir;
                    output[thin_dir] = di_thin_dir;
                    output[3]        = flux_dir;
                    return output;
                };
                
                // Boundary fluxes
                grid::cell_idx_t i_cell_bdy;
                i_cell_bdy.lb() = lb;
                i_cell_bdy.i(thin_dir)  = tile_id[thin_dir]*tile_size;
                i_cell_bdy.i(seq_dir )  = 2*tile_id[seq_dir]*tile_size;
                i_cell_bdy.i(fuse_dir)  = 2*tile_id[fuse_dir]*tile_size;
                
                const auto idx = get_offsts(inner_raw);
                const int flux_dir = idx[3];
                
                i_cell_bdy.i(0) += idx[0];
                i_cell_bdy.i(1) += idx[1];
                i_cell_bdy.i(2) += idx[2];
                
                grid::face_idx_t i_face = grid::cell_to_face(i_cell_bdy, flux_dir, 1);
                omni::retrieve(grid_img, q_img, i_face, input0);
                flux_type flux0 = flux_func(input0);
                flux0  *= inv_dx[flux_dir];
                
                
                #pragma unroll
                for (int dir = 0; dir < dim; ++dir)
                {
                    if (dir == flux_dir)
                    {
                        sh_rhs(idx[0], idx[1], idx[2]) -= flux0;
                    }
                    threads.sync();
                }
                
                #pragma unroll
                for (int t_id = 0; t_id < 2; ++t_id)
                {
                    grid::cell_idx_t i_cell;
                    i_cell.lb()         = lb;
                    i_cell.i(thin_dir)  = tile_id[thin_dir ]*tile_size        + inner_raw[thin_dir];
                    i_cell.i(seq_dir )  = (2*tile_id[seq_dir]+t_id)*tile_size + inner_raw[seq_dir];
                    i_cell.i(fuse_dir)  = 2*tile_id[fuse_dir]*tile_size       + inner_raw[fuse_dir];
                    
                    
                    grid::face_idx_t i_face;
                    
                    #pragma unroll
                    for (int idir = 0; idir < dim; ++idir)
                    {
                        i_face = grid::cell_to_face(i_cell, idir, 0);
                        omni::retrieve(grid_img, q_img, i_face, input);
                        
                        flux_type flux0 = flux_func(input);
                        flux0  *= inv_dx[idir];
                        
                        auto ii = inner_raw;
                        ii[seq_dir] += tile_size*t_id;
                        auto iil = ii;
                        iil[idir]--;
                        threads.sync();
                        sh_rhs(ii[0], ii[1], ii[2]) += flux0;
                        threads.sync();
                        if (iil[idir] >= 0) sh_rhs(iil[0], iil[1], iil[2]) -= flux0;
                    }
                }
                
                #pragma unroll
                for (int t_id = 0; t_id < 2; ++t_id)
                {
                    grid::cell_idx_t i_cell;
                    i_cell.lb()         = lb;
                    i_cell.i(thin_dir)  = tile_id[thin_dir ]*tile_size        + inner_raw[thin_dir];
                    i_cell.i(seq_dir )  = (2*tile_id[seq_dir]+t_id)*tile_size + inner_raw[seq_dir];
                    i_cell.i(fuse_dir)  = 2*tile_id[fuse_dir]*tile_size       + inner_raw[fuse_dir];
                    auto ii = inner_raw;
                    ii[seq_dir] += tile_size*t_id;
                    threads.sync();
                    flux_type& my_rhs = sh_rhs(ii[0], ii[1], ii[2]);
                    rhs_img.set_elem(i_cell, my_rhs);
                }
            });
        };
        
        // Execute the bulk workload
        dispatch::execute(outer_range, loop, kpool, k_shmem);
    }
}