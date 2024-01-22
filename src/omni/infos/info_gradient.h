#pragma once

#include "omni/info.h"

namespace spade::omni
{
    namespace info
    {
        struct gradient : public info_base<gradient> // floating-point types only, and need to be able to specify the order-of-accuracy later
        {

            constexpr static bool requires_direction = false;
            constexpr static bool is_shmem_buffered  = true;
            
            template <typename array_t, const grid::array_centering center>
            using array_data_type
            = ctrs::array<typename array_t::alias_type, array_t::grid_type::coord_point_type::size()>;
            
            template <typename grid_view_t, typename array_t, typename index_t>
            requires((index_t::centering_type() == grid::face_centered) && (grid::cell_centered == array_t::centering_type()))
            _sp_hybrid static void compute(
                const grid_view_t& ar_grid,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                using grid_t = utils::remove_all<decltype(ar_grid)>::type;
                const int idir0 = idx.dir();
                const ctrs::array<int,3> idir(idir0, (idir0+1)%ar_grid.dim(), (idir0+2)%ar_grid.dim());
                const ctrs::array<typename grid_t::coord_type, 3> invdx
                (
                    1.0/ar_grid.get_dx(0, idx.lb()), //todo: update with block index
                    1.0/ar_grid.get_dx(1, idx.lb()),
                    1.0/ar_grid.get_dx(2, idx.lb())
                );
                grid::cell_idx_t ic = grid::face_to_cell(idx, 0);

                out = 0.0;
                
                auto apply_coeff_at = [&](const int& iset, const typename array_t::value_type& coeff, const grid::cell_idx_t& idx_in)
                {
                    auto q = array.get_elem(idx_in);
                    if constexpr (ctrs::basic_array<typename array_t::alias_type>)
                    {
                        for (std::size_t i = 0; i < out[iset].size(); ++i) out[iset][i] += coeff*q[i]*invdx[iset];
                    }
                    else
                    {
                        out[iset] += coeff*q*invdx[iset];
                    }
                };
                
                apply_coeff_at(idir[0],  -1.0, ic);
                for (int ii = 1; ii < grid_t::dim(); ++ii)
                {
                    ic[idir[ii]] += 1;
                    apply_coeff_at(idir[ii],  0.25, ic);
                    ic[idir[ii]] -= 2;
                    apply_coeff_at(idir[ii], -0.25, ic);
                    ic[idir[ii]] += 1;
                }
                ic[idir[0]] += 1;
                apply_coeff_at(idir[0],   1.0, ic);
                for (int ii = 1; ii < grid_t::dim(); ++ii)
                {
                    ic[idir[ii]] += 1;
                    apply_coeff_at(idir[ii],  0.25, ic);
                    ic[idir[ii]] -= 2;
                    apply_coeff_at(idir[ii], -0.25, ic);
                    ic[idir[ii]] += 1;
                }
                using real_type = typename array_t::value_type;
                static_assert(std::same_as<typename utils::remove_all<decltype(ar_grid.get_coord_sys())>::type, coords::identity<real_type>>, "gradient transformation required for general coordinates");
                // const auto& x_face = ar_grid.get_comp_coords(iface);
                // transform_gradient(ar_grid, iface, out);
            }

            template <typename grid_view_t, typename array_t, typename index_t>
            requires((index_t::centering_type() == grid::cell_centered) && (grid::cell_centered == array_t::centering_type()))
            _sp_hybrid static void compute(
                const grid_view_t& ar_grid,
                const array_t& array,
                const index_t& idx,
                array_data_type<array_t, index_t::centering_type()>& out)
            {
                using grid_t = utils::remove_all<decltype(ar_grid)>::type;
                const ctrs::array<typename grid_t::coord_type, 3> invdx
                (
                    1.0/ar_grid.get_dx(0), //todo: update with block index
                    1.0/ar_grid.get_dx(1),
                    1.0/ar_grid.get_dx(2)
                );

                auto ic = idx;
                
                out = 0.0;
                for (int dir = 0; dir < ar_grid.dim(); ++dir)
                {
                    ic.i(dir) += 1;
                    out[dir]  += array.get_elem(ic);
                    ic.i(dir) -= 2;
                    out[dir]  -= array.get_elem(ic);
                    ic.i(dir) += 1;
                    out[dir]  *= 0.5*invdx[dir];
                }
                using real_type = typename array_t::value_type;
                static_assert(std::same_as<typename utils::remove_all<decltype(ar_grid.get_coord_sys())>::type, coords::identity<real_type>>, "gradient transformation required for general coordinates");
                // const auto& x_face = ar_grid.get_comp_coords(iface);
                // transform_gradient(ar_grid, iface, out);
            }
        };
    }
}