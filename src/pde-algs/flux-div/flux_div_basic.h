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
    static void flux_div_basic(
        const sol_arr_t& prims,
        rhs_arr_t& rhs,
        const flux_func_t& flux_func,
        const traits_t& traits)
    {
        
        // Note: this is a naive implementation for the gpu
        // at the moment, the CPU one above is faster by 2x for CPU
        using real_type        = sol_arr_t::value_type;
        using alias_type       = sol_arr_t::alias_type;
        using omni_union_type  = omni::combine_omni_stencils<flux_func_t>;
        using flux_out_t       = rhs_arr_t::alias_type;
        
        using namespace sym::literals;
        const auto& incr = algs::get_trait(traits, "pde_increment"_sym, increment);
        using incr_mode_t = typename utils::remove_all<decltype(incr)>::type;
        constexpr bool is_incr_mode = incr_mode_t::increment_mode;
        if constexpr (!is_incr_mode) rhs = real_type(0.0);
        
        
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
}