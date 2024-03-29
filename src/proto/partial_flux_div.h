#pragma once

#include <concepts>

#include "grid/grid.h"
#include "omni/omni.h"

namespace spade::proto
{    
    template <
        grid::multiblock_array sol_arr_t,
        grid::multiblock_array rhs_arr_t,
        typename flux_func_t>
    requires
        grid::has_centering_type<sol_arr_t, grid::cell_centered>
    static void flux_div_partial(
        const device::shared_vector<std::size_t>& flux_blocks,
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
        
        var_range.lower.lb() = 0UL;
        var_range.upper.lb() = flux_blocks.size();
        
        const auto fbimg = utils::make_vec_image(flux_blocks.data(prims.device()));
        
        auto load = _sp_lambda (const grid::cell_idx_t& icell_in) mutable
        {
            grid::cell_idx_t icell = icell_in;
            icell.lb() = fbimg[icell.lb()];
            
            const auto xyz = geom_image.get_comp_coords(icell);
            const real_type jac = coords::calc_jacobian(geom_image.get_coord_sys(), xyz, icell);
            auto relem = rhs_img.get_elem(icell);
            algs::static_for<0,grid_dim>([&](const auto& iidir)
            {
                constexpr int idir     = iidir.value;
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
            });
            rhs_img.set_elem(icell, relem);
        };
        dispatch::execute(var_range, load);
    }
    
}