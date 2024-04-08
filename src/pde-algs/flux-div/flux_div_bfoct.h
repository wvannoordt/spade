#pragma once

#include <concepts>

#include "grid/grid.h"
#include "omni/omni.h"

namespace spade::pde_algs
{
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t,
        typename traits_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    inline void flux_div_bfoct(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func,
        const traits_t&)
    {
        using real_type     = sol_arr_t::value_type;
        using alias_type    = sol_arr_t::alias_type;
        using omni_type     = typename flux_func_t::omni_type;
        using flux_type     = rhs_arr_t::alias_type;
        using grid_type     = typename sol_arr_t::grid_type;
        
        static_assert(std::same_as<typename grid_type::coord_sys_type, coords::identity<typename grid_type::coord_type>>, "flux divergence does not yet account for the jacobian!");
        
        const auto& grid    = prims.get_grid();
        const auto grid_img = grid.image(partition::local, prims.device());
        const auto q_img    = prims.image();
        auto rhs_img        = rhs.image();
        
        constexpr int dim = grid.dim();
        
        auto nx     = grid.get_num_cells();
        auto ngs    = prims.get_num_exchange();
        auto ntiles = nx;
        
        
        constexpr int left_extent   = -static_math::moddiv<omni_type::template min_extent<0>, 2>::value;
        constexpr int right_extent  =  static_math::moddiv<omni_type::template max_extent<0>, 2>::value + 1;
        constexpr int num_sten_vals = left_extent + right_extent;
        constexpr int ng            =  2;
        static_assert((left_extent <= ng) && (right_extent <= ng), "flux_div doesn't yet work for stencils wider than 2");
        constexpr int tile_size     =  4;
        constexpr int half_tile     = tile_size/2;
        static_assert(2*half_tile == tile_size, "Tile size must be an even number");
        
        int total_tiles = 1;
        for (int d = 0; d < nx.size(); ++d)
        {
            if ((nx[d] % tile_size) != 0)
            {
                throw except::sp_exception("flux_div requires a block size that is a multiple of 4");
            }
            if (ngs[d] != 2)
            {
                throw except::sp_exception("flux_div requires all exchange cells are size 2");
            }
            ntiles[d]    = utils::i_div_up(nx[d], tile_size);
            total_tiles *= ntiles[d];
        }
        
        using input_type = omni::stencil_data_t<omni_type, sol_arr_t>;
            
        const auto tile_range = dispatch::ranges::make_range(0, tile_size, 0, tile_size, 0, tile_size);
        dispatch::kernel_threads_t kpool(tile_range, prims.device());
        using threads_type = decltype(kpool);
            
        // NOTE: at some point we ought to add in the ability to buffer
        // RHS data to shmem for the increment operation.
        std::size_t total_sh_vals    = tile_size*tile_size*tile_size;
            
        auto k_shmem     = dispatch::shmem::make_shmem(dispatch::shmem::vec<alias_type>(total_sh_vals));
        using shmem_type = decltype(k_shmem);
            
        const auto outer_range = dispatch::ranges::make_range(0, ntiles[0]*ntiles[1]*ntiles[2], 0, int(grid.get_num_local_blocks()));
        auto loop = [=] _sp_hybrid (const ctrs::array<int, 2>& outer_raw, const threads_type& threads, shmem_type& shmem) mutable
        {
            int tile_id_1d = outer_raw[0];
            auto& shmem_vec = shmem[0_c];
            // This is potentially slow!
            // ix + nx*iy + nx*ny*iz
            ctrs::array<int, 3> tile_id;
            tile_id[0]  = tile_id_1d % ntiles[0];
            tile_id_1d -= tile_id[0];
            tile_id_1d /= ntiles[0];
            // tile_id_1d  = i0; //Does not work
            tile_id[1]  = tile_id_1d % ntiles[1];
            tile_id_1d -= tile_id[1];
            tile_id_1d /= ntiles[1];
            // tile_id_1d  = i1; //Does not work
            tile_id[2]  = tile_id_1d;
            ctrs::array<input_type, dim> input;
            int lb = outer_raw[1];
            const auto inv_dx = grid_img.get_inv_dx(lb);
            constexpr bool has_gradient = omni_type::template info_at<omni::offset_t<0,0,0>>::template contains<omni::info::gradient>;
            constexpr bool has_face_val = omni_type::template info_at<omni::offset_t<0,0,0>>::template contains<omni::info::value>;
            auto tile_data = utils::make_vec_image(shmem_vec, tile_size, tile_size, tile_size);
            threads.exec([&](const ctrs::array<int, 3>& inner_raw)
            {
                if constexpr (has_gradient)
                {
                    for (int idir = 0; idir < input.size(); ++idir)
                    {
                        //Initialize gradient to zero
                        auto& grad = omni::access<omni::info::gradient>(input[idir].root());
                        grad = real_type(0.0);
                    }
                }
                
                grid::cell_idx_t i_cell(tile_id[0]*tile_size + inner_raw[0], tile_id[1]*tile_size + inner_raw[1], tile_id[2]*tile_size + inner_raw[2], lb);
                
                #pragma unroll
                for (int quad = 0; quad < 8; ++quad)
                {
                    ctrs::array<int, 3> quad_offst{((quad & 1) > 0)?1:-1, ((quad & 2) > 0)?1:-1, ((quad & 4) > 0)?1:-1};
                    auto i_buff = i_cell;
                    
                    #pragma unroll
                    for (int d = 0; d < 3; ++d) i_buff.i(d) += quad_offst[d]*half_tile;
                    
                    // Load the data from the current adjacent quadrant
                    tile_data(inner_raw[0], inner_raw[1], inner_raw[2]) = q_img.get_elem(i_buff);
                    threads.sync();
                    
                    #pragma unroll
                    for (int idir = 0; idir < dim; ++idir)
                    {
                        
                        //Calculate partial contribution to the gradient from this quadrant
                        if constexpr (has_gradient)
                        {
                            auto& grad = omni::access<omni::info::gradient>(input[idir].root());
                            
                            // Compute all three direction components of the flux
                            int idir0 = idir + 1;
                            int idir1 = idir + 2;
                            if (idir0 >= dim) idir0 -= dim;
                            if (idir1 >= dim) idir1 -= dim;
                            
                            const auto& cond_apply = [&](int d_idir, int d_idir0, int d_idir1, const real_type coeff, int comp)
                            {
                                ctrs::array<int, 3> app_idx = inner_raw;
                                app_idx[idir]  += d_idir;
                                app_idx[idir0] += d_idir0;
                                app_idx[idir1] += d_idir1;
                                
                                app_idx[0] -= quad_offst[0]*half_tile;
                                app_idx[1] -= quad_offst[1]*half_tile;
                                app_idx[2] -= quad_offst[2]*half_tile;
                                
                                alias_type val = real_type(0.0);
                                bool loaded    = true;
                                loaded = loaded && (app_idx[0] >= 0);
                                loaded = loaded && (app_idx[0] < tile_size);
                                loaded = loaded && (app_idx[1] >= 0);
                                loaded = loaded && (app_idx[1] < tile_size);
                                loaded = loaded && (app_idx[2] >= 0);
                                loaded = loaded && (app_idx[2] < tile_size);
                                if (loaded)
                                {
                                    val = tile_data(app_idx[0], app_idx[1], app_idx[2]);
                                }
                                
                                grad[comp] += coeff*val;
                            };
                            
                            #pragma unroll
                            for (int pm = 0; pm < 2; ++pm)
                            {
                                // pm = 0 --> lower
                                // pm = 1 --> upper
                                int sig    = 2*pm - 1;
                                int d_idir = pm - 1;
                                cond_apply(d_idir,  0,  0, sig*real_type( 1.0),  idir);
                                cond_apply(d_idir, -1,  0,     real_type(-0.25), idir0);
                                cond_apply(d_idir,  1,  0,     real_type( 0.25), idir0);
                                cond_apply(d_idir,  0, -1,     real_type(-0.25), idir1);
                                cond_apply(d_idir,  0,  1,     real_type( 0.25), idir1);
                            }
                        }
                        
                        // Move the stencil values onto the stencil data from those that are loaded currently
                        auto& input_loc = input[idir];
                        ctrs::array<int, 3> ld_idx = inner_raw;
                        ld_idx[idir] -= left_extent;
                        
                        ld_idx[0] -= quad_offst[0]*half_tile;
                        ld_idx[1] -= quad_offst[1]*half_tile;
                        ld_idx[2] -= quad_offst[2]*half_tile;
                        
                        algs::static_for<0, num_sten_vals>([&](const auto& jj)
                        {
                            constexpr int j            = jj.value;
                            constexpr int offset_num   = 2*(j - left_extent) + 1;
                            using offst_t              = omni::offset_t<offset_num, 0, 0>;
                            constexpr int cell_idx     = omni::index_of<omni_type, offst_t>;
                            auto& cell_val             = omni::access<omni::info::value>(input_loc.cell(udci::idx_const_t<cell_idx>()));
                            
                            bool loaded    = true;
                            loaded = loaded && (ld_idx[0] >= 0);
                            loaded = loaded && (ld_idx[0] < tile_size);
                            loaded = loaded && (ld_idx[1] >= 0);
                            loaded = loaded && (ld_idx[1] < tile_size);
                            loaded = loaded && (ld_idx[2] >= 0);
                            loaded = loaded && (ld_idx[2] < tile_size);
                            if (loaded)
                            {
                                cell_val = tile_data(ld_idx[0], ld_idx[1], ld_idx[2]);
                            }
                            ld_idx[idir]++;
                        });
                    }
                    threads.sync();
                }
                
                if constexpr (has_gradient)
                {
                    #pragma unroll
                    for (int idir = 0; idir < dim; ++idir)
                    {
                        auto& grad = omni::access<omni::info::gradient>(input[idir].root());
                        for (int dd = 0; dd < dim; ++dd) grad[dd] *= inv_dx[dd];
                    }
                }
                
                ctrs::array<flux_type, dim> all3;
                
                #pragma unroll
                for (int idir = 0; idir < dim; ++idir)
                {
                    auto& input_loc = input[idir];
                    if constexpr (has_face_val)
                    {
                        auto& qface = omni::access<omni::info::value>(input_loc.root());
                        qface = real_type(0.0);
                        using offst0_t              = omni::offset_t<-1, 0, 0>;
                        constexpr int cell_idx0     = omni::index_of<omni_type, offst0_t>;
                        
                        using offst1_t              = omni::offset_t<1, 0, 0>;
                        constexpr int cell_idx1     = omni::index_of<omni_type, offst1_t>;
                        
                        const auto& q0 = omni::access<omni::info::value>(input_loc.cell(udci::idx_const_t<cell_idx0>()));
                        const auto& q1 = omni::access<omni::info::value>(input_loc.cell(udci::idx_const_t<cell_idx1>()));
                        qface += q0;
                        qface += q1;
                        qface *= real_type(0.5);
                    }
                    
                    const auto excluded = omni::info_list_t<omni::info::value, omni::info::gradient>();
                    const auto i_face = grid::cell_to_face(i_cell, idir, 0);
                    omni::retrieve(grid_img, q_img, i_face, input_loc, excluded);
                    
                    all3[idir] = flux_func(input_loc);
                    all3[idir] *= inv_dx[idir];
                }
                
                // Old "strategy"
                // for (int idir = 0; idir < dim; ++idir)
                // {
                //     auto i_cellL = i_cell;
                //     i_cellL.i(idir)--;
                //     if (i_cellL.i(idir) >= 0)
                //     {
                //         rhs_img.decr_elem(i_cellL, all3[idir]);
                //     }
                //     threads.sync();
                //     rhs_img.incr_elem(i_cell, all3[idir]);
                //     threads.sync();
                // }
                
                //New strategy: shmem
                auto rhs_data = utils::vec_img_cast<flux_type>(tile_data);
                
                threads.sync();
                auto& self_flx = rhs_data(inner_raw[0], inner_raw[1], inner_raw[2]);
                self_flx = rhs_img.get_elem(i_cell);
                
                threads.sync();
                #pragma unroll
                for (int idir = 0; idir < dim; ++idir)
                {
                    bool valid  = !(inner_raw[idir] == 0);
                    auto idd   = inner_raw;
                    idd[idir] -= int(valid);
                    auto& left_flx = rhs_data(idd[0], idd[1], idd[2]);
                    if (valid) left_flx -= all3[idir];
                    threads.sync();
                }
                
                #pragma unroll
                for (int idir = 0; idir < dim; ++idir)
                {
                    self_flx += all3[idir];
                }
                
                rhs_img.set_elem(i_cell, rhs_data(inner_raw[0], inner_raw[1], inner_raw[2]));
                threads.sync();
            });
        };
        
        dispatch::execute(outer_range, loop, kpool, k_shmem);
        
        
        // block boundary cleanup
        int combine_dim = -1;
        for (int d = 0; d < dim; ++d)
        {
            if (ntiles[d]%2 == 0) { combine_dim = d; break; }
        }
        if (combine_dim < 0)
        {
            throw except::sp_exception("flux_div currently requires that at least one grid dimension is a multiple of 8");
        }
        
        ctrs::array<int, 3> orange_dims = ntiles;
        orange_dims[combine_dim] /= 2;
        
        ctrs::array<int, 3> irange_dims = tile_size;
        irange_dims[combine_dim] *= 2;
        
        const auto i_range   = dispatch::ranges::make_range(0, tile_size, 0, tile_size, 0, 2);
        dispatch::kernel_threads_t kpool_cor(i_range, prims.device());
        const auto outer_range_cor = dispatch::ranges::make_range(0, orange_dims[0]*orange_dims[1]*orange_dims[2], 0, int(grid.get_num_local_blocks()));
        using threads_type_cor  = decltype(kpool_cor);
        
        using data_type = omni::stencil_data_t<omni_type, sol_arr_t>;
        
        auto loop_cor = [=] _sp_hybrid (const ctrs::array<int, 2>& outer_raw, const threads_type_cor& threads) mutable
        {
            int tile_id_1d = outer_raw[0];
            ctrs::array<int, 3> btile_id;
            ctrs::array<int, 3> ntiles_loc = ntiles;
            ntiles_loc[combine_dim] /= 2;
            btile_id[0]  = tile_id_1d % ntiles_loc[0];
            tile_id_1d -= btile_id[0];
            tile_id_1d /= ntiles_loc[0];
            // tile_id_1d  = i0; //Does not work
            btile_id[1]  = tile_id_1d % ntiles_loc[1];
            tile_id_1d -= btile_id[1];
            tile_id_1d /= ntiles_loc[1];
            // tile_id_1d  = i1; //Does not work
            btile_id[2]  = tile_id_1d;
            int lb = outer_raw[1];
            threads.exec([&](const ctrs::array<int, 3>& is_i)
            {
                auto tile_id = btile_id;
                tile_id[combine_dim] *= 2;
                tile_id[combine_dim] += is_i[2];
                
                #pragma unroll
                for (int idir = 0; idir < dim; ++idir)
                {
                    int idir0 = idir + 1;
                    int idir1 = idir + 2;
                    if (idir0 >= dim) idir0 -= dim;
                    if (idir1 >= dim) idir1 -= dim;
                    
                    grid::cell_idx_t  upper;
                    upper.lb()     =  lb;
                    upper.i()      =  tile_id[0]*tile_size;
                    upper.j()      =  tile_id[1]*tile_size;
                    upper.k()      =  tile_id[2]*tile_size;
                    upper.i(idir)  += tile_size - 1;
                    upper.i(idir0) += is_i[0];
                    upper.i(idir1) += is_i[1];
                    
                    const auto uface = grid::cell_to_face(upper, idir, 1);
                    
                    data_type input;
                    omni::retrieve(grid_img, q_img, uface, input);
                    flux_type flux = flux_func(input);
                    const auto inv_dx = grid_img.get_inv_dx(idir, upper.lb());
                    flux *= inv_dx;
                    rhs_img.decr_elem(upper, flux);
                    threads.sync();
                }
            });
        };
        dispatch::execute(outer_range_cor, loop_cor, kpool_cor);
    }
}