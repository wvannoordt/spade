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
    inline void flux_div_longf(
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
        constexpr int tspow        =  2;
        constexpr int tile_size    =  1 << tspow;
        constexpr int tsmask       = tile_size - 1;
        
        // Compute number of tiles in each direction
        int total_tiles = 1;

        
        // Create data structure for the flux kernel
        using input_type = typename omni::stencil_data_t<omni_type, sol_arr_t>;
        
        constexpr int bspow0 = 3;
        constexpr int bspow1 = 3;
        constexpr int bspow2 = 2;
        
        constexpr int bspow_big0 = bspow0 - tspow;
        constexpr int bspow_big1 = bspow1 - tspow;
        constexpr int bspow_big2 = bspow2 - tspow;
        
        constexpr int tsize0 = 1 << bspow_big0;
        constexpr int tsize1 = 1 << bspow_big1;
        constexpr int tsize2 = 1 << bspow_big2;
        
        constexpr int tmask0 = tsize0 - 1;
        constexpr int tmask1 = tsize1 - 1;
        constexpr int tmask2 = tsize2 - 1;
        
        constexpr int bsize0 = 1 << bspow0;
        constexpr int bsize1 = 1 << bspow1;
        constexpr int bsize2 = 1 << bspow2;
        
        constexpr int bmask0 = bsize0 - 1;
        constexpr int bmask1 = bsize1 - 1;
        constexpr int bmask2 = bsize2 - 1;
        
        const spade::ctrs::array<int, 3> irange_dims{bsize0, bsize1, bsize2};
        
        for (int d = 0; d < nx.size(); ++d)
        {
            ntiles[d]    = utils::i_div_up(nx[d], irange_dims[d]);
            total_tiles *= ntiles[d];
        }
        
        // Create a range for the individual tile size
        const auto tile_range = dispatch::ranges::make_range(0, irange_dims[0]*irange_dims[1]*irange_dims[2]);
        
        // Create collaborative threads object over the tile range
        dispatch::kernel_threads_t kpool(tile_range, prims.device());
        using threads_type = decltype(kpool);
        
        const std::size_t total_sh_vals = 0;//tile_size*tile_size*tile_size;
        
        // Create shared memory array
        auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<flux_type>(total_sh_vals));
        using shmem_type = decltype(k_shmem);
        
        // Create a grid-range over the number of tiles in each block and the number of blocks
        const auto outer_range = dispatch::ranges::make_range(0, int(grid.get_num_local_blocks()));
        // Create a lambda for the bulk workload
        
        auto loop = [=] _sp_hybrid (const int& outer_raw, const threads_type& threads, shmem_type& shmem) mutable
        {
            // Compute the tile index for this thread block
            // This is the input data structure for the flux calculation
            input_type input;
            
            // Block index
            int lb = outer_raw;
            
            // 1 / (grid spacing) for differentiation
            const auto inv_dx_native = grid_img.get_inv_dx(lb);
            
            // Convert 1/deltaX into the necessary precision for this computation
            ctrs::array<real_type, dim> inv_dx;
            #pragma unroll
            for (int d = 0; d < dim; ++d) inv_dx[d] = inv_dx_native[d];
            
            threads.exec([&](const int& inner_raw)
            {
                for (int k_blk = 0; k_blk < ntiles[2]; ++k_blk)
                {
                    for (int j_blk = 0; j_blk < ntiles[1]; ++j_blk)
                    {
                        for (int i_blk = 0; i_blk < ntiles[0]; ++i_blk)
                        {
                            
                            grid::cell_idx_t i_cell;
                            
                            int idx = inner_raw;
                            int k_small = idx & tsmask;
                            idx = idx >> tspow;
                            int j_small = idx & tsmask;
                            idx = idx >> tspow;
                            int i_small = idx & tsmask;
                            idx = idx >> tspow;
                            
                            int k_big = idx & tmask2;
                            idx = idx >> bspow_big2;
                            int j_big = idx & tmask1;
                            idx = idx >> bspow_big1;
                            int i_big = idx & tmask0;
                            idx = idx >> bspow_big0;
                            
                            i_cell.lb() = lb;
                            i_cell.i()  = irange_dims[0]*i_blk + i_big*tile_size + i_small;
                            i_cell.j()  = irange_dims[1]*j_blk + j_big*tile_size + j_small;
                            i_cell.k()  = irange_dims[2]*k_blk + k_big*tile_size + k_small;
                            
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
                            
                            grid::face_idx_t i_face;
                            algs::static_for<0,dim>([&](const auto& ii)
                            {
                                constexpr int idir = ii.value;
                                algs::static_for<0,2>([&](const auto& jj)
                                {
                                    constexpr int t_pm = jj.value;
                                    i_face = grid::cell_to_face(i_cell, idir, t_pm);
                                    omni::retrieve(grid_img, q_img, i_face, input);
                                    
                                    flux_type flux0 = flux_func(input);
                                    flux0 *= inv_dx[idir];
                                    const auto coeff = 1 - 2*t_pm;
                                    my_rhs += coeff*flux0;
                                });
                            });
                            
                            rhs_img.set_elem(i_cell, my_rhs);
                        }
                    }
                }
            });
        };
        
        // Execute the bulk workload
        dispatch::execute(outer_range, loop, kpool, k_shmem);
    }
}
    
