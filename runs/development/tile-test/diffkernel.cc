#include "diffkernel.h"

void diffkernel(int i, int j, int k, int lb, array_t& rhs, const array_t& q)
{
    const auto& grid = rhs.get_grid();
    auto dx = grid.get_block_box(lb).size(0)/grid.get_num_cells(0);
    auto dy = grid.get_block_box(lb).size(1)/grid.get_num_cells(1);
    auto dz = grid.get_block_box(lb).size(2)/grid.get_num_cells(2);
    const double a0 = -1.0/12.0;
    const double a1 =  4.0/3.0;
    const double a2 = -5.0/2.0;
    const double a3 =  4.0/3.0;
    const double a4 = -1.0/12.0;
    rhs(i, j, k, lb) += (a0*q(i+2,   j,   k, lb) + a1*q(i+1,   j,   k, lb) + a2*q(i, j, k, lb) + a3*q(i-1,   j,   k, lb) + a4*q(i-2,   j,   k, lb))/(dx*dx);
    rhs(i, j, k, lb) += (a0*q(i,   j+2,   k, lb) + a1*q(i,   j+1,   k, lb) + a2*q(i, j, k, lb) + a3*q(i,   j-1,   k, lb) + a4*q(i,   j-2,   k, lb))/(dy*dy);
    rhs(i, j, k, lb) += (a0*q(i,   j,   k+2, lb) + a1*q(i,   j,   k+1, lb) + a2*q(i, j, k, lb) + a3*q(i,   j,   k-1, lb) + a4*q(i,   j,   k-2, lb))/(dz*dz);
}