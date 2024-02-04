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
            
            spade::dispatch::kernel_threads_t kpool(dispatch::ranges::make_range(0, num_faces), prims.device());
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
                
                auto input_data = omni::make_buffered_data<omni_type, sol_arr_t>(shmem, -1000);
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
}