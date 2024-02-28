#pragma once

#include <concepts>

#include "grid/grid.h"
#include "omni/omni.h"

namespace spade::pde_algs
{
    template <typename flux_func_t> concept is_flux_functor =
    requires(const flux_func_t& f)
    {
        //TODO
        f;
    };
    
    template <typename arr_t, typename flux_eval_t> concept has_flux_compatibility =
    requires(const arr_t& t0, const flux_eval_t& t1)
    {
        //TODO
        t0;
    };

    enum block_flux_policy
    {
        block_flux_all,
        block_flux_interior
    };
    
    /*
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const bound_box_t<bool,sol_arr_t::grid_type::dim()>& domain_boundary_flux,
        const flux_funcs_t&... flux_funcs)
    {
        using real_type  = sol_arr_t::value_type;
        using alias_type = sol_arr_t::alias_type;
        const grid::multiblock_grid auto& ar_grid = prims.get_grid();
        const auto geom_image = ar_grid.image(prims.device());
        
        const auto prims_img = prims.image();
        auto rhs_img         = rhs.image();

        auto block_range = range(0,ar_grid.get_num_local_blocks());

        using omni_union_type = omni::combine_omni_stencils<flux_funcs_t...>;

        for (auto lb_idir: block_range*range(0, ar_grid.dim()))
        {
            int  lb_loc  = lb_idir[0];
            auto lb_glob = ar_grid.get_partition().to_global(utils::tag[partition::local](lb_loc));
            int idir = lb_idir[1];
            bound_box_t<int,3> flux_bounds;
            for (auto d:range(0,3))
            {
                flux_bounds.min(d) = (ar_grid.dim()==2 && d==2)?0:-1;
                flux_bounds.max(d) = ar_grid.get_num_cells(d);
            }
            const auto& iboundary = ar_grid.is_domain_boundary(lb_glob);
            if (iboundary.min(idir) && !domain_boundary_flux.min(idir)) flux_bounds.min(idir) = 0;
            if (iboundary.max(idir) && !domain_boundary_flux.max(idir)) flux_bounds.max(idir) = ar_grid.get_num_cells(idir) - 1;
            auto g0 = range(flux_bounds.min(0),flux_bounds.max(0));
            auto g1 = range(flux_bounds.min(1),flux_bounds.max(1));
            auto g2 = range(flux_bounds.min(2),flux_bounds.max(2));

            auto grid_range = g0*g1*g2;
            for (auto idx: grid_range)
            {
                auto lb_loc_des = utils::tag[partition::local](lb_loc);
                const real_type dx = geom_image.get_dx(idir, lb_loc);
                grid::cell_idx_t il(idx[0_c], idx[1_c], idx[2_c], lb_loc);
                grid::cell_idx_t ir(idx[0_c], idx[1_c], idx[2_c], lb_loc);
                ir[idir] += 1;
                grid::face_idx_t iface = grid::cell_to_face(il, idir, 1);
                const auto xyz_comp_l = geom_image.get_comp_coords(il);
                const auto xyz_comp_r = geom_image.get_comp_coords(ir);
                const real_type jac_l = coords::calc_jacobian(geom_image.get_coord_sys(), xyz_comp_l, il);
                const real_type jac_r = coords::calc_jacobian(geom_image.get_coord_sys(), xyz_comp_r, ir);
                
                using omni_type = omni_union_type;
                using data_type = omni::stencil_data_t<omni_type, sol_arr_t>;
                data_type data;
                omni::retrieve(geom_image, prims_img, iface, data);
                using flux_out_t = rhs_arr_t::alias_type;
                flux_out_t accum(0.0);
                utils::foreach_param([&](const auto& flux_func)
                {   
                    using flux_func_t = decltype(flux_func);
                    using omni_type   = utils::remove_all<flux_func_t>::type::omni_type;

                    auto flux_data = omni::interpret_stencil<omni_type>(data);
                    accum += flux_func(flux_data);
                    
                }, flux_funcs...);

                //todo: fix this garbage
                if constexpr (ctrs::basic_array<alias_type>)
                {
                    for (int n = 0; n < accum.size(); ++n)
                    {
                        rhs_img(n, il) -= jac_l*accum[n]/(dx);
                        rhs_img(n, ir) += jac_r*accum[n]/(dx);
                    }
                }
                else
                {
                    rhs_img(il) -= jac_l*accum/(dx);
                    rhs_img(ir) += jac_r*accum/(dx);
                }
            }
        }
    }
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename... flux_funcs_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_funcs_t&... flux_funcs)
    {
        const std::size_t dim = sol_arr_t::grid_type::dim();
        bound_box_t<bool,dim> domain_boundary_flux;
        for (auto i: range(0,dim))
        {
            domain_boundary_flux.max(i) = true;
            domain_boundary_flux.min(i) = true;
        }
        flux_div(prims, rhs, domain_boundary_flux, flux_funcs...);
    }
    */

    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div_OLD(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func)
    {
        
        // Note: this is a naive implementation for the gpu
        // at the moment, the CPU one above is faster by 2x for CPU
        using real_type        = sol_arr_t::value_type;
        using alias_type       = sol_arr_t::alias_type;
        using omni_union_type  = omni::combine_omni_stencils<flux_func_t>;
        using flux_out_t       = rhs_arr_t::alias_type;
        
        const auto& ar_grid    = prims.get_grid();
        const auto geom_image  = ar_grid.image(prims.device());
        const auto prims_img   = prims.image();
        auto rhs_img           = rhs.image();
        
        constexpr int grid_dim = ar_grid.dim();
        auto var_range         = dispatch::support_of(prims, grid::exclude_exchanges);
        
        for(int idir = 0; idir < ar_grid.dim(); ++idir)
        {
            auto load = [=] _sp_hybrid (const grid::cell_idx_t& icell) mutable
            {
                const auto xyz = geom_image.get_comp_coords(icell);
                const real_type jac = coords::calc_jacobian(geom_image.get_coord_sys(), xyz, icell);
                auto relem = rhs_img.get_elem(icell);
                
                const real_type inv_dx = real_type(1.0)/geom_image.get_dx(idir, icell.lb());
                algs::static_for<0,2>([&](const auto& ipm)
                {
                    constexpr int pm = ipm.value;
                    auto iface = grid::cell_to_face(icell, idir, pm);
                    
                    using omni_type = omni_union_type;
                    using data_type = omni::stencil_data_t<omni_type, sol_arr_t>;
                    
                    data_type data;
                    omni::retrieve(geom_image, prims_img, iface, data);
                    
                    flux_out_t accum = 0.0;
                    auto flux_data = omni::interpret_stencil<omni_type>(data);
                    accum += flux_func(flux_data);
                    accum *= jac*(real_type(1.0)-real_type(2.0)*pm)*inv_dx;
                    relem += accum;
                });
                rhs_img.set_elem(icell, relem);
            };
            dispatch::execute(var_range, load);
        }
    }
  
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func)
    {
        
        using real_type        = sol_arr_t::value_type;
        using alias_type       = sol_arr_t::alias_type;
        using omni_type        = typename flux_func_t::omni_type;
        using flux_out_t       = rhs_arr_t::alias_type;
        
        const auto& ar_grid    = prims.get_grid();
        const auto geom_image  = ar_grid.image(partition::local, prims.device());
        const auto prims_img   = prims.image();
        auto rhs_img           = rhs.image();
        
        constexpr int grid_dim = ar_grid.dim();
        
        for(int idir = 0; idir < ar_grid.dim(); ++idir)
        {
            int idir0 = (idir+1)%grid_dim;
            int idir1 = (idir+2)%grid_dim;
            
            const int num_faces = ar_grid.get_num_cells(idir) + 1;
            
            auto k_shmem = omni::get_shmem<omni_type>(num_faces, prims);
            using shmem_type = decltype(k_shmem);
            
            dispatch::kernel_threads_t kpool(dispatch::ranges::make_range(0, num_faces), prims.device());
            auto range = dispatch::ranges::make_range(0, ar_grid.get_num_cells(idir0), 0, ar_grid.get_num_cells(idir1), 0, ar_grid.get_num_local_blocks());
            using index_type = decltype(range)::index_type;
            using threads_type = decltype(kpool);
            
            auto loop = [=] _sp_hybrid (const index_type& idx, const threads_type& threads, shmem_type& shmem) mutable
            {
                // d, i, j, k, lb
                grid::face_idx_t iface;
                iface.dir()    = idir;
                iface.i(idir0) = idx[0];
                iface.i(idir1) = idx[1];
                iface.lb()     = idx[2];
                
                auto input_data = omni::make_buffered_data<omni_type, sol_arr_t>(shmem, -1000, 0, 1, 0, 1);
                threads.exec([&](const int i_face_line)
                {
                    //Need to come up with something more elegant here
                    iface.i(idir)       = i_face_line;
                    input_data.line_idx = i_face_line;
                    omni::retrieve_shared(geom_image, prims_img, iface, input_data, shmem);
                });
                
                threads.sync();
                
                threads.exec([&](const int i_face_line)
                {
                    iface.i(idir)       = i_face_line;
                    input_data.line_idx = i_face_line;
                    const real_type inv_dx = real_type(1.0)/geom_image.get_dx(idir, iface.lb());
                    
                    omni::retrieve(geom_image, prims_img, iface, input_data.supplemental);
                    // This is private, so any residual cell requiring this must be addressed
                    using output_type = typename flux_func_t::output_type;
                    output_type flux = flux_func(input_data);
                    flux *= (inv_dx);
                    
                    //Here we decrement the right-hand element
                    const auto icell_r    = face_to_cell(iface, 1);
                    const auto xyz_r      = geom_image.get_coords(icell_r);
                    const real_type jac_r = coords::calc_jacobian(geom_image.get_coord_sys(), xyz_r, icell_r);
                    if (i_face_line < num_faces - 1) rhs_img.incr_elem(icell_r, jac_r*flux);
                    
                    threads.sync();
                    
                    // Here we increment the left-hand element
                    const auto icell_l    = face_to_cell(iface, 0);
                    const auto xyz_l      = geom_image.get_coords(icell_r);
                    const real_type jac_l = coords::calc_jacobian(geom_image.get_coord_sys(), xyz_l, icell_l);
                    if (i_face_line > 0) rhs_img.decr_elem(icell_l, jac_l*flux);
                });
                
            };
            spade::dispatch::execute(range, loop, kpool, k_shmem);
        }
    }
    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    inline void flux_div_NEW(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func)
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
        
        for (int parity = 0; parity < 4; ++parity)
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
                
                if (tpar == parity)
                {
                    input_type input;
                    int lb = outer_raw[1];
                    const auto inv_dx = grid_img.get_inv_dx(lb);
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
                            if (is_interior)           rhs_img.incr_elem(i_cell,   flux);
                            threads.sync();
                            if (do_lft && is_interior) rhs_img.decr_elem(i_cell_l, flux);
                            threads.sync();
                        }
                    });
                }
            };
            dispatch::execute(outer_range, loop, kpool, k_shmem);
        }
        
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
                    const auto inv_dx = grid_img.get_inv_dx(idir, upper.lb());
                    flux *= inv_dx;
                    bool valid = (upper.i(idir0) < nx[idir0]) && (upper.i(idir1) < nx[idir1]);
                    if (valid) rhs_img.decr_elem(upper, flux);
                });
            };
            dispatch::execute(outer_range, loop, kpool);
        }
    }
}