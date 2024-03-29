#pragma once

#include <concepts>
#include "core/bounding_box.h"
#include "core/static_for.h"

namespace spade::boundary
{
    using identifier_t = bound_box_t<bool, 3>;
    static const identifier_t xmin(true,  false, false, false, false, false);
    static const identifier_t xmax(false, true,  false, false, false, false);
    static const identifier_t ymin(false, false, true,  false, false, false);
    static const identifier_t ymax(false, false, false, true,  false, false);
    static const identifier_t zmin(false, false, false, false, true,  false);
    static const identifier_t zmax(false, false, false, false, false, true );
    
    template <const int extrap_order>
    struct extrap_t
    {
        constexpr static bool is_extrap_type_v = true;
        constexpr static int order = extrap_order;
    };
    
    template <typename T> concept is_extrap_type = T::is_extrap_type_v;
    
    template <const int ord> static constexpr extrap_t<ord> extrapolate = extrap_t<ord>();
}

namespace spade::algs
{
    template <typename arr_t, typename kernel_t>
    static void boundary_fill(arr_t& arr, const boundary::identifier_t& boundaries, const kernel_t kern)
    {
        using arr_val_t = typename arr_t::value_type;
        using crd_val_t = typename arr_t::grid_type::coord_type;
        static_assert(arr_t::centering_type() == grid::cell_centered, "boundary filling only implemented for cell-centered arrays");
        const auto& grid   = arr.get_grid();
        const auto g_img   = grid.image(arr.device());
        auto arr_img = arr.image();
        for (int ibndy = 0; ibndy < 6; ++ibndy)
        {
            int idir = ibndy / 2;
            int pm   = ibndy % 2;
            if (boundaries(idir, pm))
            {
                grid::cell_idx_t ll, ur;
                ll.lb() = 0;
                ur.lb() = g_img.boundary_blocks[ibndy].size();
                for (int i = 0; i < 3; ++i)
                {
                    ll[i] = -arr.get_num_exchange(i);
                    ur[i] =  grid.get_num_cells(i) + arr.get_num_exchange(i);
                }
                
                // Lower boundary
                if (pm == 0) ur[idir] = 0;
                
                // Upper boundary
                if (pm == 1) ll[idir] = grid.get_num_cells(idir);
                
                const dispatch::ranges::grid_idx_range_t bd_range(ll, ur, arr.device());
                using index_type = decltype(bd_range)::index_type;
                
                // TODO: fix this horrific mess when NVIDIA actually implement proper device lambdas
                if constexpr (boundary::is_extrap_type<kernel_t>)
                {
                    auto extrap_bdy = _sp_lambda(const index_type& ii) mutable
                    {
                        // NOTE: ii.lb() is incorrect, need to update it with the value from the domain boundary table
                        // this is weird but necessary
                        grid::cell_idx_t fill_idx = ii;
                        fill_idx.lb() = g_img.boundary_blocks[ibndy][ii.lb()];
                        
                        using coeff_t = typename arr_t::fundamental_type;
                        using alias_t = typename arr_t::alias_type;
                        
                        alias_t ghost_value = coeff_t(0.0);
                        int i0 = pm*(g_img.num_cell[idir] - kernel_t::order - 1);
                        
                        algs::static_for<0, kernel_t::order + 1>([&](const auto& jj)
                        {
                            constexpr int j = jj.value;
                            coeff_t lj = coeff_t(1.0);
                            algs::static_for<0, kernel_t::order + 1>([&](const auto& ii)
                            {
                                constexpr int i = ii.value;
                                if constexpr (i != j)
                                {
                                    const auto x  = fill_idx[idir];
                                    const auto xj = i0 + j;
                                    const auto xi = i0 + i;
                                    lj = lj*(x - xi)/(xj - xi);
                                }
                            });
                            auto idx_interp = fill_idx;
                            idx_interp[idir] = i0 + j;
                            
                            ghost_value = ghost_value+lj*arr_img.get_elem(idx_interp);
                        });
                        arr_img.set_elem(fill_idx, ghost_value);
                    };
                    dispatch::execute(bd_range, extrap_bdy);
                }
                else
                {
                    auto fill_bdy = _sp_lambda(const index_type& ii) mutable
                    {
                        // NOTE: ii.lb() is incorrect, need to update it with the value from the domain boundary table
                        // this is weird but necessary
                        grid::cell_idx_t fill_idx = ii;
                        fill_idx.lb() = g_img.boundary_blocks[ibndy][ii.lb()];
                        grid::cell_idx_t image_idx = fill_idx;
                        if (pm == 0) image_idx[idir] = -1 - fill_idx[idir];
                        if (pm == 1) image_idx[idir] = 2*g_img.num_cell[idir] - (fill_idx[idir] + 1);
                        //auto domain_val = arr_img.get_elem(image_idx);
                        //auto ghost_val  = kern(domain_val, idir);
                        auto domain_val = arr_img.get_elem(image_idx);
                        auto fill_val   = arr_img.get_elem(fill_idx);
                        auto x_g        = g_img.get_coords(fill_idx);
                        auto kern2 = kern;
                        auto ghost_val  = [&]()
                        {
                          if constexpr (std::invocable<kernel_t, typename arr_t::alias_type, typename arr_t::alias_type, typename coords::point_t<crd_val_t>,  int>)
                            return kern2(domain_val, fill_val, x_g, idir);
                          else return kern2(domain_val, idir);
                        }();
                        arr_img.set_elem(fill_idx, ghost_val);
                    };
                    dispatch::execute(bd_range, fill_bdy);
                }
            }
        }
    }
}
