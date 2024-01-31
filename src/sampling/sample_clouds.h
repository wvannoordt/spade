#pragma once

#include "core/ctrs.h"

namespace spade::sampling
{
    static struct multilinear_t
    {
        constexpr static int stencil_size() { return 8; }
        
        template <
            ctrs::basic_array                              indices_t,
            ctrs::basic_array                              coeffs_t,
            grid::multiblock_grid                          grid_t,
            ctrs::basic_array                              reduced_t,
            ctrs::basic_array                              delta_t,
            std::invocable<typename indices_t::value_type> exclude_t
            >
        requires (
            (indices_t::size() == coeffs_t::size()) &&
            indices_t::size() >= stencil_size()
            )
        bool try_make_cloud(
            indices_t& indices,
            coeffs_t& coeffs,
            const grid_t& grid,
            const typename indices_t::value_type& landed_cell,
            const typename grid_t::coord_point_type& x_sample,
            const reduced_t& reduced_idx,
            const delta_t& deltai,
            const exclude_t& exclude_crit) const
        {
            using coeff_t = typename coeffs_t::value_type;
            using pnt_t   = typename grid_t::coord_point_type;
            int nmx = 4;
            if (grid_t::dim()==3) nmx = 8;
            
            grid::cell_idx_t lower_left = landed_cell;
            for (int d = 0; d < 3; ++d) lower_left.i(d) += (deltai[d]<0)?deltai[d]:0;
            
            indices = landed_cell;
            coeffs  = 0.0;
            
            for (int ct = 0; ct < nmx; ++ct)
            {
                ctrs::array<int, 3> di = 0;
                di[0] = ct%2;
                di[1] = (ct/2)%2;
                di[2] = (ct/4)%2;
                
                indices[ct].i()  = lower_left.i() + di[0];
                indices[ct].j()  = lower_left.j() + di[1];
                indices[ct].k()  = lower_left.k() + di[2];
                indices[ct].lb() = lower_left.lb();
                
                //Need to check these!
                coeff_t coeff = 1.0;
                for (int d = 0; d < grid_t::dim(); ++d)
                {
                    coeff_t dir_coeff0 = reduced_idx[d] - lower_left[d];
                    coeff_t dir_coeff1 = 1.0 - dir_coeff0;
                    coeff_t dir_coeff  = di[d]*dir_coeff0 + (1.0-di[d])*dir_coeff1;
                    coeff *= dir_coeff;
                }
                coeffs[ct] = coeff;
            }
            return true;
        }
    } multilinear;
}