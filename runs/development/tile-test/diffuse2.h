#pragma once

#include "diffkernel.h"

void ijk_loop(const array_t& q, array_t& rhs, int ni, int nj, int nk, int nlb, void (*fcn)(int,int,int,int,array_t&,const array_t&))
{
    for (int lb = 0; lb < nlb; lb++)
    {
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    fcn(i, j, k, lb, rhs, q);
                }
            }
        }
    }
}

static void diffuse2(array_t& rhs, const array_t& q)
{
    const auto& grid = rhs.get_grid();
    int nlb = grid.get_num_local_blocks();
    int nk = grid.get_num_cells(2);
    int nj = grid.get_num_cells(1);
    int ni = grid.get_num_cells(0);
    ijk_loop(q, rhs, ni, nj, nk, nlb, diffkernel);
}