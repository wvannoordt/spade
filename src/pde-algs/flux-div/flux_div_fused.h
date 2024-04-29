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
    inline void flux_div_fused(
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
        
        // Compute number of tiles in each direction
        int total_tiles = 1;
        for (int d = 0; d < nx.size(); ++d)
        {
            ntiles[d]    = utils::i_div_up(nx[d], tile_size);
            total_tiles *= ntiles[d];
        }
        
        // Create data structure for the flux kernel
        
        using omni_left_type  = omni::shift_stencil<omni_type, omni::offset_t< 0, 0, 0>>;
        using omni_right_type = omni::shift_stencil<omni_type, omni::offset_t< 2, 0, 0>>;
        using omni_all_type   = omni::stencil_union<omni_left_type, omni_right_type>;
        
        using omni_left_view  = omni::shift_stencil<omni_all_type, omni::offset_t< 0, 0, 0>>;
        using omni_right_view = omni::shift_stencil<omni_all_type, omni::offset_t<-2, 0, 0>>;
        
        
        using input_type     = typename omni::stencil_data_t<omni_type, sol_arr_t>;
        using input_all_type = typename omni::stencil_data_t<omni_all_type, sol_arr_t>;
        
        // debug::print_type(omni_left_view());
        // debug::print_type(omni_right_view());
        
        input_all_type inp;
        grid::face_idx_t i_facee(0, 5, 5, 5, 5);
        
        omni::access<omni::info::value>(inp.cell(0_c)) = 0.0;
        omni::access<omni::info::value>(inp.cell(1_c)) = 1.0;
        omni::access<omni::info::value>(inp.cell(2_c)) = 2.0;
        omni::access<omni::info::value>(inp.cell(3_c)) = 3.0;
        omni::access<omni::info::value>(inp.cell(4_c)) = 4.0;
        
        // omni::retrieve(grid_img, q_img, i_facee, inp);
        const auto d0L = omni::interpret_stencil_at<omni_type, omni::offset_t<0, 0, 0>>(inp);
        const auto d0R = omni::interpret_stencil_at<omni_type, omni::offset_t<2, 0, 0>>(inp);
        
        
        // print(omni::access<omni::info::value>(d0L.cell(0_c)));
        // print(omni::access<omni::info::value>(d0L.cell(1_c)));
        // print(omni::access<omni::info::value>(d0L.cell(2_c)));
        // print(omni::access<omni::info::value>(d0L.cell(3_c)));
        // print("aaaa");
        // print(omni::access<omni::info::value>(d0R.cell(0_c)));
        // print(omni::access<omni::info::value>(d0R.cell(1_c)));
        // print(omni::access<omni::info::value>(d0R.cell(2_c)));
        // print(omni::access<omni::info::value>(d0R.cell(3_c)));
        
        // debug::print_type(all_1_t());
        // debug::print_type(all_2_t());
        
        // debug::print_type(omni_left_view());
        // debug::print_type(omni_right_view());
        
        // std::cin.get();
        
        constexpr int thin_dir = 0;
        constexpr int fuse_dir = 2;
        constexpr int seq_dir  = 1;
        
        spade::ctrs::array<int, 3> irange_dims = tile_size;
        irange_dims[fuse_dir] *= 2;
        // irange_dims[seq_dir] *= 2;
        
        ntiles[fuse_dir] /= 2;
        ntiles[seq_dir]  /= 2;
        
        constexpr int n_tid_thin = 1;
        ntiles[thin_dir]  /= n_tid_thin;
        
        // Create a range for the individual tile size
        const auto tile_range = dispatch::ranges::make_range(0, irange_dims[0], 0, irange_dims[1], 0, irange_dims[2]);
        
        // Create collaborative threads object over the tile range
        dispatch::kernel_threads_t kpool(tile_range, prims.device());
        using threads_type = decltype(kpool);
        
        const std::size_t total_sh_vals = 0;//tile_size*tile_size*tile_size;
        
        // Create shared memory array
        auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<flux_type>(total_sh_vals));
        using shmem_type = decltype(k_shmem);
        
        // Create a grid-range over the number of tiles in each block and the number of blocks
        constexpr int lb_fs = 1;
        const auto outer_range = dispatch::ranges::make_range(0, ntiles[0], 0, ntiles[1]*ntiles[2], 0, int(grid.get_num_local_blocks())/lb_fs);
        
        for (int lb_par = 0; lb_par < lb_fs; ++lb_par)
        {
            // Create a lambda for the bulk workload
            auto loop = [=] _sp_hybrid (const ctrs::array<int, 3>& outer_raw, const threads_type& threads, shmem_type& shmem) mutable
            {
                // Compute the tile index for this thread block
                int tile_id_1d = outer_raw[1];
                auto& shmem_vec = shmem[0_c];
                ctrs::array<int, 3> tile_id;
                tile_id[0]  = outer_raw[0];
                tile_id[1]  = tile_id_1d % ntiles[1];
                tile_id_1d -= tile_id[1];
                tile_id_1d /= ntiles[1];
                tile_id[2]  = tile_id_1d;
                
                // Block index
                int lb = lb_fs*outer_raw[2] + lb_par;
                
                // 1 / (grid spacing) for differentiation
                const auto inv_dx_native = grid_img.get_inv_dx(lb);
                
                // Convert 1/deltaX into the necessary precision for this computation
                ctrs::array<real_type, dim> inv_dx;
                #pragma unroll
                for (int d = 0; d < dim; ++d) inv_dx[d] = inv_dx_native[d];
                
                threads.exec([&](const ctrs::array<int, 3>& inner_raw)
                {
                    constexpr int n_tid      = 2;
                    
                    // #pragma unroll << intentionally left out!
                    for (int t_id_thin = 0; t_id_thin < n_tid_thin; ++t_id_thin)
                    {
                        #pragma unroll
                        for (int t_id = 0; t_id < n_tid; ++t_id)
                        {
                            const int t_id_eff = t_id*t_id_thin + (1-t_id_thin)*(1-t_id);
                            grid::cell_idx_t i_cell;
                            i_cell.lb()         = lb;
                            i_cell.i(thin_dir)  = (n_tid_thin*tile_id[thin_dir] + t_id_thin)*tile_size        + inner_raw[thin_dir];
                            i_cell.i(seq_dir )  = (2*tile_id[seq_dir]+t_id_eff)*tile_size + inner_raw[seq_dir];
                            i_cell.i(fuse_dir)  = 2*tile_id[fuse_dir]*tile_size       + inner_raw[fuse_dir];
                            flux_type my_rhs;
                            const auto nothing = rhs_img.size();
                            constexpr bool is_incr_modetmp = is_incr_mode;
                            if constexpr (is_incr_modetmp)
                            {
                                my_rhs = rhs_img.get_elem(i_cell);
                            }
                            else
                            {
                                my_rhs = real_type(0.0);
                            }
                            
                            constexpr bool use_fused_data = true;
                            const auto ff = flux_func;
                            const auto qq = q_img;
                            if constexpr (!use_fused_data)
                            {
                                // Overwrite the RHS (residual) element
                                grid::face_idx_t i_face;
                                
                                #pragma unroll
                                for (int idir = 0; idir < dim; ++idir)
                                // algs::static_for<0, dim>([&](const auto& ii)
                                {
                                    // constexpr int idir = ii.value;
                                    #pragma unroll
                                    for (int t_pm = 0; t_pm < 2; ++t_pm)
                                    // algs::static_for<0, 2>([&](const auto& jj)
                                    {
                                        // constexpr int t_pm = jj.value;
                                        i_face = grid::cell_to_face(i_cell, idir, t_pm);
                                        input_type input;
                                        omni::retrieve(grid_img, q_img, i_face, input);
                                        
                                        flux_type flux0 = flux_func(input);
                                        flux0 *= inv_dx[idir];
                                        const auto coeff = 1 - 2*t_pm;
                                        my_rhs += coeff*flux0;
                                    // });
                                    }
                                }
                                // });
                            }
                            else
                            {
                                // #pragma unroll
                                // for (int idir = 0; idir < dim; ++idir)
                                algs::static_for<0, dim>([&](const auto& ii)
                                {
                                    constexpr int idir = ii.value;
                                    const grid::face_idx_t i_face = grid::cell_to_face(i_cell, idir, 0);
                                    input_all_type input;
                                    omni::retrieve(grid_img, q_img, i_face, input);
                                    
                                    const auto ll = omni::interpret_stencil_at<omni_type, omni::offset_t<0, 0, 0>>(input);
                                    const auto rr = omni::interpret_stencil_at<omni_type, omni::offset_t<2, 0, 0>>(input);
                                    
                                    flux_type fluxL = flux_func(ll);
                                    flux_type fluxR = flux_func(rr);
                                    
                                    fluxL -= fluxR;
                                    fluxL *= inv_dx[idir];
                                    my_rhs += fluxL;
                                });
                                // }
                            }
                            
                            rhs_img.set_elem(i_cell, my_rhs);
                        }
                    }
                });
            };
            
            // Execute the bulk workload
            dispatch::execute(outer_range, loop, kpool, k_shmem);
        }
    }
}