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
    inline void flux_div_fldbc(
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
        using input_type = typename omni::stencil_data_t<omni_type, sol_arr_t>;
        
        // Create a range for the individual tile size
        const auto tile_range = dispatch::ranges::make_range(0, tile_size, 0, tile_size, 0, tile_size);
        
        // Create collaborative threads object over the tile range
        dispatch::kernel_threads_t kpool(tile_range, prims.device());
        using threads_type = decltype(kpool);
        
        // This is the dimensions for one of the views of the shared memory array
        auto view0 = ctrs::make_array(tile_size+2*ng, tile_size, tile_size);
        
        // Get shared memory size in elements
        std::size_t num_sh_vals_flx  = view0[0]*view0[1]*view0[2];
        std::size_t total_sh_vals    = num_sh_vals_flx;
        
        // Create chared memory array
        auto k_shmem = dispatch::shmem::make_shmem(dispatch::shmem::vec<alias_type>(total_sh_vals));
        using shmem_type = decltype(k_shmem);
        
        // Create a grid-range over the number of tiles in each block and the number of blocks
        const auto outer_range = dispatch::ranges::make_range(0, ntiles[0], 0, ntiles[1]*ntiles[2], 0, int(grid.get_num_local_blocks()));
        
        // Create a lambda for the bulk workload
        auto loop = [=] _sp_hybrid (const ctrs::array<int, 3>& outer_raw, const threads_type& threads, shmem_type& shmem) mutable
        {
            // Compute the tile index for this thread block
            int tile_id_1d = outer_raw[1];
            auto& shmem_vec = shmem[0_c];
            ctrs::array<int, 3> tile_id;
            tile_id[0] = outer_raw[0];
            tile_id[1]  = tile_id_1d % ntiles[1];
            tile_id_1d -= tile_id[1];
            tile_id_1d /= ntiles[1];
            tile_id[2]  = tile_id_1d;

            // This is the input data structure for the flux calculation
            input_type input;
            
            // Block index
            int lb = outer_raw[2];
            
            // 1 / (grid spacing) for differentiation
            const auto inv_dx_native = grid_img.get_inv_dx(lb);
            
            // Convert 1/deltaX into the necessary precision for this computation
            ctrs::array<real_type, dim> inv_dx;
            #pragma unroll
            for (int d = 0; d < dim; ++d) inv_dx[d] = inv_dx_native[d];
            
            // Compile-time check if the stencil has a gradient at the face
            constexpr bool has_gradient = omni_type::template info_at<omni::offset_t<0,0,0>>::template contains<omni::info::gradient>;
            
            // Compile-time check if the stencil has a value at the face (requires averaging from cells)
            constexpr bool has_face_val = omni_type::template info_at<omni::offset_t<0,0,0>>::template contains<omni::info::value>;
            
            threads.exec([&](const ctrs::array<int, 3>& inner_raw)
            {
                // Calculate the cell index based on tile number and thread id
                grid::cell_idx_t i_cell;
                i_cell.lb()     = lb;
                i_cell.i()      = tile_id[0]*tile_size + inner_raw[0];
                i_cell.j()      = tile_id[1]*tile_size + inner_raw[1];
                i_cell.k()      = tile_id[2]*tile_size + inner_raw[2];
                
                // Mask out threads that are operating out of bounds
                bool is_interior = true;
                is_interior = is_interior && (i_cell.i() >= 0);
                is_interior = is_interior && (i_cell.j() >= 0);
                is_interior = is_interior && (i_cell.k() >= 0);
                is_interior = is_interior && (i_cell.i() < nx[0]);
                is_interior = is_interior && (i_cell.j() < nx[1]);
                is_interior = is_interior && (i_cell.k() < nx[2]);
                
                // "alias_type" is the type that is stored over the entire state array
                // this value is re-used multiple times, so is just loaded into
                // a local variable for the duration of the kernel
                alias_type my_elem;                        
                if (is_interior) my_elem = q_img.get_elem(i_cell);
                
                // "flux_type" is the type that is stored over the entire residual array
                // Every thread will buffer its own element of the residual array
                flux_type my_rhs;
                
                // Quirk for first-capture variables, nothing happens here
                const auto nothing = rhs_img.size();
                constexpr bool is_incr_modetmp = is_incr_mode;
                
                // Overwrite the RHS (residual) element
                if constexpr (is_incr_modetmp)
                {
                    if (is_interior) my_rhs = rhs_img.get_elem(i_cell);
                }
                else
                {
                    my_rhs = real_type(0.0);
                }
                
                // Direction (x,y,z) loop
                #pragma unroll
                for (int idir = 0; idir < dim; ++idir)
                {
                    // Compute the tangential directions
                    int idir0 = idir + 1;
                    int idir1 = idir + 2;
                    if (idir0 >= dim) idir0 -= dim;
                    if (idir1 >= dim) idir1 -= dim;
                    
                    // Array of tangential directions
                    auto other_dirs = ctrs::make_array(idir0, idir1);
                    
                    // face-index of the face whose flux this thread calculates
                    auto i_face    = grid::cell_to_face(i_cell, idir, 0);
                    
                    if constexpr (has_gradient)
                    {
                        // Here, we compute tangential components of the gradient
                        auto& gradient  = omni::access<omni::info::gradient>(input.root());
                        gradient = real_type(0.0);
                        
                        // Create an appropriately-sized view of the shared memory
                        auto faces_view = utils::make_vec_image(
                            shmem_vec,
                            2,         // -/+
                            2,         // tan_dir = 0 or 1
                            tile_size, // idir[norm_dir]
                            tile_size  // idir[tan_dir]
                        );
                        
                        // Compute indices for face data buffering
                        int pm        = inner_raw[idir0] & 1;
                        int fdir_id   = (inner_raw[idir0] & 2) >> 1;
                        int i_nrm     = inner_raw[idir];
                        int i_tan     = inner_raw[idir1];
                        int fdir      = other_dirs[fdir_id];
                        int other_dir = other_dirs[1-fdir_id];
                        
                        // Lower-left indices of this tile
                        auto ll = ctrs::make_array(tile_id[0]*tile_size, tile_id[1]*tile_size, tile_id[2]*tile_size);
                        
                        // Compute the cell address to load into the face array
                        auto i_targ     = i_cell;
                        i_targ.i(idir)  = i_cell.i(idir);
                        i_targ.i(idir0) = ll[idir0];
                        i_targ.i(idir1) = ll[idir1];
                        
                        i_targ.i(fdir)      += (pm*tile_size + (1-pm)*(-1));
                        i_targ.i(other_dir) += i_tan;
                        
                        // Compute the store address
                        auto& target = faces_view(pm, fdir_id, i_nrm, i_tan);
                        target = q_img.get_elem(i_targ);
                        
                        // Populate the bulk data with thread's own element
                        auto voldata = utils::make_vec_image(shmem_vec, tile_size, tile_size, tile_size);
                        voldata.ptr  = faces_view.end();
                        
                        voldata(inner_raw[0], inner_raw[1], inner_raw[2]) = my_elem;
                        threads.sync();
                        
                        // Computing the gradient in tangential directions using face and volume data
                        #pragma unroll
                        for (int i_norm_pm = 0; i_norm_pm < 2; ++i_norm_pm)
                        {
                            int idx_nrm = utils::max(i_nrm - (1 - i_norm_pm), 0);
                            bool lower_face_nrm = (i_nrm == 0);
                            
                            // Tangential direction loop
                            #pragma unroll
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
                                
                                // Differentiate in tangential direction as though all data is available
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
                                
                                // If data is coming from the edge, substitute it in
                                if (block_left)  left_val  = edge_val;
                                if (block_right) right_val = edge_val;
                                
                                gradient[dir_here] += coeff*right_val;
                                gradient[dir_here] -= coeff*left_val;
                            }
                        }
                        // Note that by this time, not all tangential gradients are computed.
                        // Cells on the lower face of the tile (in normal direction)
                        // still require a ring of data around the edge, which
                        // will be handled later
                        threads.sync();
                    } // End fringe gradient pt. 1
                    
                    // Buffering the values along the flux stencil
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
                    
                    // Check to see if the data being buffered is redundant, in
                    // which case we assign an element on the ring to be buffered
                    if (inner_raw[idir] == tile_size - 1)
                    {
                        // Calculate the index of the ring element
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
                        // Finish tangential gradient calculation using the ring information
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
                            
                            // Similar strategy as above - calculate as if the ring value is used,
                            // block out otherwise
                            auto val_left   = vals(ii_l[0], ii_l[1], ii_l[2]);
                            auto val_right  = vals(ii_r[0], ii_r[1], ii_r[2]);
                            
                            if (block_left)  val_left  = edge_val;
                            if (block_right) val_right = edge_val;
                            
                            auto coeff          = real_type(0.25)*inv_dx[real_dir]*(inner_raw[idir]==0);
                            gradient[real_dir] += coeff*val_right;
                            gradient[real_dir] -= coeff*val_left;
                        }
                    }
                    
                    
                    // Now we fill the stencil values from the volume data that is loaded
                    constexpr int num_stencil_vals = 2*ng;
                    auto ii_l = inner_raw;
                    
                    // Loop over stencil elements
                    algs::static_for<0, num_stencil_vals>([&](const auto& iidx)
                    {
                        const auto idxxx = udci::idx_const_t<iidx.value>();
                        auto& st_val = omni::access<omni::info::value>(input.cell(idxxx));
                        st_val = vals(ii_l[0], ii_l[1], ii_l[2]);
                        ii_l[idir]++;
                    });
                    
                    threads.sync();
                    
                    // Compute face-average value directly from stencil
                    if constexpr (has_face_val)
                    {
                        auto& face_val    = omni::access<omni::info::value>(input.root());
                        constexpr int lft = omni::index_of<omni_type, omni::offset_t<-1, 0, 0>>;
                        constexpr int rgt = omni::index_of<omni_type, omni::offset_t< 1, 0, 0>>;
                        auto& q_upper     = omni::access<omni::info::value>(input.cell(udci::idx_const_t<rgt>()));
                        auto& q_lower     = omni::access<omni::info::value>(input.cell(udci::idx_const_t<lft>()));
                        face_val          = q_upper;
                        face_val         += q_lower;
                        face_val         *= real_type(0.5);
                    }
                    
                    // Compute normal gradient directly from stencil
                    if constexpr (has_gradient)
                    {
                        auto& gradient  = omni::access<omni::info::gradient>(input.root());
                        constexpr int lft = omni::index_of<omni_type, omni::offset_t<-1, 0, 0>>;
                        constexpr int rgt = omni::index_of<omni_type, omni::offset_t< 1, 0, 0>>;
                        auto& q_upper     = omni::access<omni::info::value>(input.cell(udci::idx_const_t<rgt>()));
                        auto& q_lower     = omni::access<omni::info::value>(input.cell(udci::idx_const_t<lft>()));
                        gradient[idir]    = q_upper;
                        gradient[idir]   -= q_lower;
                        gradient[idir]   *= inv_dx[idir];
                    }
                    
                    // Check if stencil has any other necessary information to load, in which case
                    // simply do it using the naive approach
                    const auto excluded = omni::info_list_t<omni::info::value, omni::info::gradient>();
                    if (is_interior) omni::retrieve(grid_img, q_img, i_face, input, excluded);
                    
                    // Calculate the flux (expensive compute workload)
                    flux_type flux = flux_func(input);
                    flux *= inv_dx[idir];

                    auto i_cell_l = i_cell;
                    i_cell_l.i(idir)--;
                    
                    // Don't modify left residual cell if at the tile boundary
                    bool do_lft = i_cell_l.i(idir) >= 0;
                    do_lft = do_lft && (inner_raw[idir] > 0);
                    
                    
                    // Create a view of the shared memory vector to store the residual elements
                    auto tmp     = utils::make_vec_image(shmem_vec, tile_size, tile_size, tile_size);
                    auto rawdata = utils::vec_img_cast<flux_type>(tmp);
                    auto i_rhs_mod = inner_raw;
                    
                    // Add the flux to thread's own residual
                    my_rhs += flux;
                    
                    // Load this threads RHS into shared memory
                    rawdata(i_rhs_mod[0], i_rhs_mod[1], i_rhs_mod[2]) = my_rhs;
                    threads.sync();
                    
                    // Subtract the flux from the neighbor's RHS element
                    i_rhs_mod[idir]--;                            
                    if (do_lft && is_interior) rawdata(i_rhs_mod[0], i_rhs_mod[1], i_rhs_mod[2]) -= flux;
                    threads.sync();
                    
                    // Pull the new RHS element from shared memory
                    my_rhs = rawdata(inner_raw[0], inner_raw[1], inner_raw[2]);
                    threads.sync();
                } // Direction loop
                
                if (is_interior) rhs_img.set_elem(i_cell, my_rhs);
                
            });
        };
        
        // Execute the bulk workload
        dispatch::execute(outer_range, loop, kpool, k_shmem);

        // Block boundary cleanup
        // We compute two tiles at once, compute which direction to "fuse" blocks in
        int combine_dim = -1;
        for (int d = 0; d < dim; ++d)
        {
            if (ntiles[d]%2 == 0) { combine_dim = d; break; }
        }
        if (combine_dim < 0)
        {
            throw except::sp_exception("flux_div currently requires that at least one grid dimension is a multiple of 8");
        }
        
        // Computing block / grid dimensions
        ctrs::array<int, 3> orange_dims = ntiles;
        orange_dims[combine_dim] /= 2;
        
        ctrs::array<int, 3> irange_dims = tile_size;
        irange_dims[combine_dim] *= 2;
        
        // Compute block range
        const auto i_range   = dispatch::ranges::make_range(0, tile_size, 0, tile_size, 0, 2);
        dispatch::kernel_threads_t kpool_cor(i_range, prims.device());
        
        // Compute grid range
        const auto outer_range_cor = dispatch::ranges::make_range(0, orange_dims[0]*orange_dims[1]*orange_dims[2], 0, int(grid.get_num_local_blocks()));
        using threads_type_cor  = decltype(kpool_cor);
        
        using data_type = omni::stencil_data_t<omni_type, sol_arr_t>;
            
        int nblocks = grid.get_num_local_blocks();
        ctrs::array<int, 3> ntiles_loc = ntiles;
        ntiles_loc[combine_dim] /= 2;
        
        // Flux correction workload
        auto loop_cor = [=] _sp_hybrid (const ctrs::array<int, 2>& outer_raw, const threads_type_cor& threads) mutable
        {
            int tile_id_1d = outer_raw[0];
            ctrs::array<int, 3> btile_id;
            btile_id[0]  = tile_id_1d % ntiles_loc[0];
            tile_id_1d -= btile_id[0];
            tile_id_1d /= ntiles_loc[0];
            // tile_id_1d  = i0; //Does not work
            btile_id[1]  = tile_id_1d % ntiles_loc[1];
            tile_id_1d -= btile_id[1];
            tile_id_1d /= ntiles_loc[1];
            // tile_id_1d  = i1; //Does not work
            btile_id[2]  = tile_id_1d;
            // int lb = outer_raw[1];
            int lb = nblocks - 1 - outer_raw[1];
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
                    
                    // Naive flux implementation for correction
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