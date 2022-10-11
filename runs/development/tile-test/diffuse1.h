#pragma once
static void diffuse1(auto& rhs, const auto& q)
{
    const auto& grid = q.get_grid();
    int nlb = grid.get_num_local_blocks();
    int nk = grid.get_num_cells(2);
    int nj = grid.get_num_cells(1);
    int ni = grid.get_num_cells(0);
    const real_t a0 = -1.0/12.0;
    const real_t a1 =  4.0/3.0;
    const real_t a2 = -5.0/2.0;
    const real_t a3 =  4.0/3.0;
    const real_t a4 = -1.0/12.0;
    for (int lb = 0; lb < nlb; ++lb)
    {
        auto dx = grid.get_block_box(lb).size(0)/ni;
        auto dy = grid.get_block_box(lb).size(1)/nj;
        auto dz = grid.get_block_box(lb).size(2)/nk;
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    rhs(i, j, k, lb) += (a0*q(i+2,   j,   k, lb) + a1*q(i+1,   j,   k, lb) + a2*q(i, j, k, lb) + a3*q(i-1,   j,   k, lb) + a4*q(i-2,   j,   k, lb))/(dx*dx);
                    rhs(i, j, k, lb) += (a0*q(i,   j+2,   k, lb) + a1*q(i,   j+1,   k, lb) + a2*q(i, j, k, lb) + a3*q(i,   j-1,   k, lb) + a4*q(i,   j-2,   k, lb))/(dy*dy);
                    rhs(i, j, k, lb) += (a0*q(i,   j,   k+2, lb) + a1*q(i,   j,   k+1, lb) + a2*q(i, j, k, lb) + a3*q(i,   j,   k-1, lb) + a4*q(i,   j,   k-2, lb))/(dz*dz);
                }
            }
        }
    }
}