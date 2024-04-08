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
        const bool use_parity,
        typename traits_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    inline void flux_div_ldbal(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func,
        const tldbal_t<use_parity>&,
        const traits_t& traits)
    {
        using real_type     = sol_arr_t::value_type;
        using alias_type    = sol_arr_t::alias_type;
        using omni_type     = typename flux_func_t::omni_type;
        using flux_type     = rhs_arr_t::alias_type;
        using grid_type     = typename sol_arr_t::grid_type;
        
        static_assert(std::same_as<typename grid_type::coord_sys_type, coords::identity<typename grid_type::coord_type>>, "flux divergence does not yet account for the jacobian!");
        
        using namespace sym::literals;
        
        const auto& incr = algs::get_trait(traits, "pde_increment"_sym, increment);
        using incr_mode_t = typename utils::remove_all<decltype(incr)>::type;
        constexpr bool is_incr_mode = incr_mode_t::increment_mode;
        
        const auto& grid    = prims.get_grid();
        const auto grid_img = grid.image(partition::local, prims.device());
        const auto q_img    = prims.image();
        auto rhs_img        = rhs.image();
        
        constexpr int dim = grid.dim();
        
        auto nx     = grid.get_num_cells();
        auto ngs    = prims.get_num_exchange();
        auto ntiles = nx;
        
        constexpr int left_extent  = -static_math::moddiv<omni_type::template min_extent<0>, 2>::value;
        constexpr int right_extent =  static_math::moddiv<omni_type::template max_extent<0>, 2>::value + 1;
        constexpr int ng           =  utils::max(left_extent, right_extent);
        constexpr int tile_size    =  4;
        
        int total_tiles = 1;
        for (int d = 0; d < nx.size(); ++d)
        {
            ntiles[d]    = utils::i_div_up(nx[d], tile_size);
            total_tiles *= ntiles[d];
        }
        
        using input_type_no_ptr = omni::stencil_data_t<omni_type, sol_arr_t>;
        using input_type_ptr    = omni::stencil_data_t<omni_type, sol_arr_t>;
        
        constexpr bool vals_in_shmem = false;
        
        using input_type = typename utils::choose<vals_in_shmem, input_type_ptr, input_type_no_ptr>;
        
        constexpr bool use_parity_loop = use_parity;
        constexpr int nparity = use_parity_loop ? 4 : 1;
        for (int parity = 0; parity < nparity; ++parity)
        {
            const auto tile_range = dispatch::ranges::make_range(0, tile_size, 0, tile_size, 0, tile_size);
            dispatch::kernel_threads_t kpool(tile_range, prims.device());
            using threads_type = decltype(kpool);
            
            // NOTE: at some point we ought to add in the ability to buffer
            // RHS data to shmem for the increment operation.
            
            // View for the other fluxes
            auto view0 = ctrs::make_array(tile_size+2*ng, tile_size, tile_size);
            
            std::size_t num_sh_vals_flx  = view0[0]*view0[1]*view0[2];
            std::size_t total_sh_vals    = num_sh_vals_flx;
            
            auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<alias_type>(total_sh_vals));
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
                
                int p0 = (tile_id[0] & 1) + 2*(tile_id[1] & 1) + 4*(tile_id[2] & 1);
                int p1 = ~p0 & 7;
                
                int tpar = utils::min(p0, p1);
                
                if ((tpar == parity) || (!use_parity_loop))
                {
                    input_type input;
                    int lb = outer_raw[1];
                    const auto inv_dx_native = grid_img.get_inv_dx(lb);
                    ctrs::array<real_type, dim> inv_dx;
                    for (int d = 0; d < dim; ++d) inv_dx[d] = inv_dx_native[d];
                    constexpr bool has_gradient = omni_type::template info_at<omni::offset_t<0,0,0>>::template contains<omni::info::gradient>;
                    constexpr bool has_face_val = omni_type::template info_at<omni::offset_t<0,0,0>>::template contains<omni::info::value>;
                    
                    threads.exec([&](const ctrs::array<int, 3>& inner_raw)
                    {
                        grid::cell_idx_t i_cell;
                        i_cell.lb()     = lb;
                        i_cell.i()      = tile_id[0]*tile_size + inner_raw[0];
                        i_cell.j()      = tile_id[1]*tile_size + inner_raw[1];
                        i_cell.k()      = tile_id[2]*tile_size + inner_raw[2];
                        
                        bool is_interior = true;
                        is_interior = is_interior && (i_cell.i() >= 0);
                        is_interior = is_interior && (i_cell.j() >= 0);
                        is_interior = is_interior && (i_cell.k() >= 0);
                        is_interior = is_interior && (i_cell.i() < nx[0]);
                        is_interior = is_interior && (i_cell.j() < nx[1]);
                        is_interior = is_interior && (i_cell.k() < nx[2]);
                        
                        alias_type my_elem;                        
                        if (is_interior) my_elem = q_img.get_elem(i_cell);
                        
                        flux_type my_rhs;
                        if (is_interior) my_rhs = rhs_img.get_elem(i_cell);
                        
                        constexpr bool is_incr_modetmp = is_incr_mode;
                        if constexpr (!is_incr_modetmp) my_rhs = real_type(0.0);
                        
                        #pragma unroll
                        for (int idir = 0; idir < dim; ++idir)
                        {
                            int idir0 = idir + 1;
                            int idir1 = idir + 2;
                            if (idir0 >= dim) idir0 -= dim;
                            if (idir1 >= dim) idir1 -= dim;
                            
                            auto other_dirs = ctrs::make_array(idir0, idir1);
                            
                            auto i_face    = grid::cell_to_face(i_cell, idir, 0);
                            if constexpr (has_gradient)
                            {
                                //Fringe gradient calculation
                                auto& gradient  = omni::access<omni::info::gradient>(input.root());
                                gradient = real_type(0.0);
                                
                                auto faces_view = utils::make_vec_image(
                                    shmem_vec,
                                    2,         // -/+
                                    2,         // tan_dir = 0 or 1
                                    tile_size, // idir[norm_dir]
                                    tile_size  // idir[tan_dir]
                                );
                                int pm        = inner_raw[idir0] & 1;
                                int fdir_id   = (inner_raw[idir0] & 2) >> 1;
                                int i_nrm     = inner_raw[idir];
                                int i_tan     = inner_raw[idir1];
                                int fdir      = other_dirs[fdir_id];
                                int other_dir = other_dirs[1-fdir_id];
                                
                                
                                auto ll = ctrs::make_array(tile_id[0]*tile_size, tile_id[1]*tile_size, tile_id[2]*tile_size);
                                
                                // Compute the buffer address
                                auto i_targ     = i_cell;
                                i_targ.i(idir)  = i_cell.i(idir);
                                i_targ.i(idir0) = ll[idir0];
                                i_targ.i(idir1) = ll[idir1];
                                
                                i_targ.i(fdir)      += (pm*tile_size + (1-pm)*(-1));
                                i_targ.i(other_dir) += i_tan;
                                
                                auto& target = faces_view(pm, fdir_id, i_nrm, i_tan);
                                target = q_img.get_elem(i_targ);
                                threads.sync();
                                
                                auto voldata = utils::make_vec_image(shmem_vec, tile_size, tile_size, tile_size);
                                voldata.ptr  = faces_view.end();
                                
                                voldata(inner_raw[0], inner_raw[1], inner_raw[2]) = my_elem;
                                threads.sync();
                                
                                #pragma unroll
                                for (int i_norm_pm = 0; i_norm_pm < 2; ++i_norm_pm)
                                {
                                    int idx_nrm = utils::max(i_nrm - (1 - i_norm_pm), 0);
                                    bool lower_face_nrm = (i_nrm == 0);
                                    for (int i_td = 0; i_td < 2; ++i_td)
                                    {
                                        int dir_here       = other_dirs[i_td];
                                        int other_dir_here = other_dirs[1-i_td];
                                        int pm_here        = inner_raw[dir_here] > 0;
                                        int tdir_idx       = inner_raw[other_dir_here];
                                        real_type coeff    = real_type(0.25)*inv_dx[dir_here];
                                        
                                        // Zero out if the data isn't actually available
                                        if (i_norm_pm == 0 && lower_face_nrm) coeff = real_type(0.0);
                                        const auto& edge_val = faces_view(pm_here, i_td, idx_nrm, tdir_idx);
                                        auto raws = inner_raw;
                                        raws[idir] += (i_norm_pm - 1);
                                        bool block_left  = raws[dir_here] == 0;
                                        bool block_right = raws[dir_here] == (tile_size - 1);
                                        raws[dir_here]--;
                                        raws[dir_here] = utils::max(0, raws[dir_here]);
                                        auto left_val  = voldata(raws[0], raws[1], raws[2]);
                                        raws[dir_here] = inner_raw[dir_here];
                                        raws[dir_here]++;
                                        raws[dir_here] = utils::min(tile_size-1, raws[dir_here]);
                                        auto right_val = voldata(raws[0], raws[1], raws[2]);
                                        if (block_left)  left_val  = edge_val;
                                        if (block_right) right_val = edge_val;
                                        
                                        gradient[dir_here] += coeff*right_val;
                                        gradient[dir_here] -= coeff*left_val;
                                    }
                                }
                                threads.sync();
                            } // End fringe gradient pt. 1
                            
                            // Now, we buffer the flowfield values
                            auto sizes           = ctrs::make_array(tile_size, tile_size, tile_size);
                            sizes[idir]         += 2*ng - 1;
                            auto vals            = utils::make_vec_image(shmem_vec, sizes);
                            auto ring            = utils::make_vec_image(shmem_vec, 2, 2, tile_size); //pm, face_dir, idx
                            ring.ptr             = vals.end();
                            auto ii              = inner_raw;
                            ii[idir]            += ng;
                            
                            threads.sync();
                            
                            auto buf_cell = i_cell;
                            vals(ii[0], ii[1], ii[2]) = my_elem;
                            threads.sync();
                            int other_buf_offset = -(1 + 2*inner_raw[idir]);
                            if (inner_raw[idir] == 2) other_buf_offset = 2;
                            buf_cell.i(idir) += other_buf_offset;
                            ii[idir] += other_buf_offset;
                            alias_type* to_set = &vals(ii[0], ii[1], ii[2]);
                            if (inner_raw[idir] == tile_size - 1)
                            {
                                int ring_pm        = inner_raw[idir0] & 1;
                                int ring_fdir_id   = (inner_raw[idir0] & 2) >> 1;
                                int ring_fdir      = other_dirs[ring_fdir_id];
                                int ring_other_dir = other_dirs[1 - ring_fdir_id];
                                int ring_idx       = inner_raw[idir1];
                                buf_cell.i()  = tile_id[0]*tile_size;
                                buf_cell.j()  = tile_id[1]*tile_size;
                                buf_cell.k()  = tile_id[2]*tile_size;
                                buf_cell.i(idir)--;
                                buf_cell.i(ring_other_dir) += ring_idx;
                                buf_cell.i(ring_fdir)      += ring_pm*tile_size + (1-ring_pm)*(-1);
                                to_set = &ring(ring_pm, ring_fdir_id, ring_idx);
                            }
                            *to_set = q_img.get_elem(buf_cell);
                            threads.sync();
                            
                            if constexpr (has_gradient)
                            {
                                // Finish fringe calculation
                                auto& gradient  = omni::access<omni::info::gradient>(input.root());
                                for (int fdir_id = 0; fdir_id < 2; ++fdir_id)
                                {
                                    int real_dir         = other_dirs[fdir_id];
                                    int other_dir        = other_dirs[1-fdir_id];
                                    int pm_here          = int(inner_raw[real_dir]>0);
                                    const auto& edge_val = ring(pm_here, fdir_id, inner_raw[other_dir]);
                                    
                                    bool block_left  = inner_raw[real_dir]==0;
                                    bool block_right = inner_raw[real_dir]==tile_size-1;
                                    auto ii_r        = inner_raw;
                                    ii_r[idir]      += ng;
                                    ii_r[idir]--;
                                    
                                    auto ii_l = ii_r;
                                    ii_r[real_dir]++;
                                    ii_l[real_dir]--;
                                    
                                    ii_r[real_dir] = utils::min(ii_r[real_dir], tile_size - 1);
                                    ii_l[real_dir] = utils::max(ii_l[real_dir], 0);
                                    auto val_left   = vals(ii_l[0], ii_l[1], ii_l[2]);
                                    auto val_right  = vals(ii_r[0], ii_r[1], ii_r[2]);
                                    
                                    if (block_left)  val_left  = edge_val;
                                    if (block_right) val_right = edge_val;
                                    
                                    auto coeff          = real_type(0.25)*inv_dx[real_dir]*(inner_raw[idir]==0);
                                    gradient[real_dir] += coeff*val_right;
                                    gradient[real_dir] -= coeff*val_left;
                                }
                                
                                // Do the normal gradient calculation
                                auto i_upper     = inner_raw;
                                i_upper[idir]   += ng;
                                auto& q_upper    = vals(i_upper[0], i_upper[1], i_upper[2]);
                                i_upper[idir]   -= 1;
                                auto& q_lower    = vals(i_upper[0], i_upper[1], i_upper[2]);
                                gradient[idir]  += q_upper;
                                gradient[idir]  -= q_lower;
                                gradient[idir]  *= inv_dx[idir];
                            }
                            
                            if constexpr (has_face_val)
                            {
                                auto& face_val   = omni::access<omni::info::value>(input.root());
                                face_val         = real_type(0.0);
                                auto i_upper     = inner_raw;
                                i_upper[idir]   += ng;
                                auto& q_upper    = vals(i_upper[0], i_upper[1], i_upper[2]);
                                auto tmp0 = i_upper;
                                i_upper[idir]   -= 1;
                                auto tmp1 = i_upper;
                                auto& q_lower    = vals(i_upper[0], i_upper[1], i_upper[2]);
                                face_val        += q_upper;
                                face_val        += q_lower;
                                face_val        *= real_type(0.5);
                            }
                            
                            // assign the stencil values
                            constexpr int num_stencil_vals = 2*ng;
                            auto ii_l = inner_raw;
                            ii_l[idir] += ng; //Just trust me
                            ii_l[idir] -= ng;
                            algs::static_for<0, num_stencil_vals>([&](const auto& iidx)
                            {
                                const auto idxxx = udci::idx_const_t<iidx.value>();
                                auto& st_val = omni::access<omni::info::value>(input.cell(idxxx));
                                st_val = vals(ii_l[0], ii_l[1], ii_l[2]);
                                ii_l[idir]++;
                            });
                            
                            const auto excluded = omni::info_list_t<omni::info::value, omni::info::gradient>();
                            if (is_interior) omni::retrieve(grid_img, q_img, i_face, input, excluded);
                            
                            flux_type flux = flux_func(input);
                            flux *= inv_dx[idir];
                            auto i_cell_l = i_cell;
                            i_cell_l.i(idir)--;
                            
                            bool do_lft = i_cell_l.i(idir) >= 0;
                            if constexpr (!use_parity_loop) do_lft = do_lft && (inner_raw[idir] > 0);
                            threads.sync();
                            
                            //Residual modification
                            auto tmp     = utils::make_vec_image(shmem_vec, tile_size, tile_size, tile_size);
                            auto rawdata = utils::vec_img_cast<flux_type>(tmp);
                            auto i_rhs_mod = inner_raw;
                            my_rhs += flux;
                            rawdata(i_rhs_mod[0], i_rhs_mod[1], i_rhs_mod[2]) = my_rhs;
                            threads.sync();
                            i_rhs_mod[idir]--;                            
                            if (do_lft && is_interior) rawdata(i_rhs_mod[0], i_rhs_mod[1], i_rhs_mod[2]) -= flux;
                            threads.sync();
                            my_rhs = rawdata(inner_raw[0], inner_raw[1], inner_raw[2]);
                            threads.sync();
                        }
                        
                        if (is_interior) rhs_img.set_elem(i_cell, my_rhs);
                        
                    });
                }
            };
            dispatch::execute(outer_range, loop, kpool, k_shmem);
        }
        
        if constexpr (use_parity_loop)
        {
            #pragma unroll
            for (int idir = 0; idir < dim; ++idir)
            {
                constexpr int correc_tile_size = 4;
                int idir0 = idir + 1;
                int idir1 = idir + 2;
                if (idir0 >= dim) idir0 -= dim;
                if (idir1 >= dim) idir1 -= dim;
                const auto i_range = dispatch::ranges::make_range(0, correc_tile_size, 0, correc_tile_size);
                const int nt_0 = utils::i_div_up(nx[idir0], correc_tile_size);
                const int nt_1 = utils::i_div_up(nx[idir1], correc_tile_size);
                dispatch::kernel_threads_t kpool(i_range, prims.device());
                const auto outer_range = dispatch::ranges::make_range(0, nt_0, 0, nt_1, 0, int(grid.get_num_local_blocks()));
                using threads_type     = decltype(kpool);
                
                using data_type = omni::stencil_data_t<omni_type, sol_arr_t>;
                
                auto loop = [=] _sp_hybrid (const ctrs::array<int, 3>& is_o, const threads_type& threads) mutable
                {
                    threads.exec([&](const ctrs::array<int, 2>& is_i)
                    {
                        grid::cell_idx_t upper;
                        upper.lb()     = is_o[2];
                        upper.i(idir)  = nx[idir]-1;
                        upper.i(idir0) = is_o[0]*correc_tile_size + is_i[0];
                        upper.i(idir1) = is_o[1]*correc_tile_size + is_i[1];
                        
                        const auto uface = grid::cell_to_face(upper, idir, 1);
                        
                        data_type input;
                        omni::retrieve(grid_img, q_img, uface, input);
                        flux_type flux = flux_func(input);
                        const auto inv_dx = real_type(grid_img.get_inv_dx(idir, upper.lb()));
                        flux *= inv_dx;
                        bool valid = (upper.i(idir0) < nx[idir0]) && (upper.i(idir1) < nx[idir1]);
                        if (valid) rhs_img.decr_elem(upper, flux);
                    });
                };
                dispatch::execute(outer_range, loop, kpool);
            }
        }
        else
        {
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
                        const auto inv_dx = real_type(grid_img.get_inv_dx(idir, upper.lb()));
                        flux *= inv_dx;
                        rhs_img.decr_elem(upper, flux);
                        threads.sync();
                    }
                });
            };
            dispatch::execute(outer_range_cor, loop_cor, kpool_cor);
        }
    }
}