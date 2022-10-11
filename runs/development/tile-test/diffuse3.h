#pragma once


template <
    typename index_t,
    const std::size_t index_rank,
    typename fcn_t,
    const std::size_t index_val
    >
requires (index_val < 0)
static void exec_loop(
    const spade::bound_box_t<index_t, index_rank>& bds,
    const spade::ctrs::array<index_t, index_rank>& index,
    const fcn_t& fcn)
{
    fcn(index);
}

template <
    typename index_t,
    const std::size_t index_rank,
    typename fcn_t,
    const std::size_t index_val
    >
requires (index_val >= 0)
static void exec_loop(
    const spade::bound_box_t<index_t, index_rank>& bds,
    spade::ctrs::array<index_t, index_rank>& index,
    const fcn_t& fcn)
{
    auto imin = bds.min(index_val);
    auto imax = bds.max(index_val);
    for (index_t i = imin; i < imax; ++i)
    {
        index[index_val] = i;
        exec_loop<index_t, index_rank, fcn_t, index_val-1>(bds, index, fcn);
    }
}

template <
    typename index_t,
    const std::size_t index_rank,
    typename fcn_t
    >
static void multi_loop(
    const spade::bound_box_t<index_t, index_rank> bds,
    const fcn_t& fcn)
{
    spade::ctrs::array<index_t, index_rank> index;
    exec_loop<index_t, index_rank, fcn_t, index_rank-1>(bds, index, fcn);
}

static void diffuse3(auto& rhs, const auto& q)
{
    const auto& grid = q.get_grid();
    int nlb = grid.get_num_local_blocks();
    using v3i = spade::ctrs::array<int, 3>;
    spade::bound_box_t<int,3> nijk = 0;
    nijk.max(0) = grid.get_num_cells(0);
    nijk.max(1) = grid.get_num_cells(1);
    nijk.max(2) = grid.get_num_cells(2);
    const real_t a0 = -1.0/12.0;
    const real_t a1 =  4.0/3.0;
    const real_t a2 = -5.0/2.0;
    const real_t a3 =  4.0/3.0;
    const real_t a4 = -1.0/12.0;
    for (int lb = 0; lb < nlb; ++lb)
    {
        auto dx = grid.get_block_box(lb).size(0)/nijk.max(0);
        auto dy = grid.get_block_box(lb).size(1)/nijk.max(1);
        auto dz = grid.get_block_box(lb).size(2)/nijk.max(2);
        auto fcn = [&](const v3i& ijk) -> void
        {
            int i = ijk[0];
            int j = ijk[1];
            int k = ijk[2];
            rhs(i, j, k, lb) += (a0*q(i+2,   j,   k, lb) + a1*q(i+1,   j,   k, lb) + a2*q(i, j, k, lb) + a3*q(i-1,   j,   k, lb) + a4*q(i-2,   j,   k, lb))/(dx*dx);
            rhs(i, j, k, lb) += (a0*q(i,   j+2,   k, lb) + a1*q(i,   j+1,   k, lb) + a2*q(i, j, k, lb) + a3*q(i,   j-1,   k, lb) + a4*q(i,   j-2,   k, lb))/(dy*dy);
            rhs(i, j, k, lb) += (a0*q(i,   j,   k+2, lb) + a1*q(i,   j,   k+1, lb) + a2*q(i, j, k, lb) + a3*q(i,   j,   k-1, lb) + a4*q(i,   j,   k-2, lb))/(dz*dz);
        };
        multi_loop(nijk, fcn);
    }
}