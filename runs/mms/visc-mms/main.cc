#include "spade.h"

int main(int argc, char** argv)
{
    typedef double real_t;
    typedef spade::fluid_state::prim_t<real_t> prim_t;
    typedef spade::fluid_state::flux_t<real_t> flux_t;
    typedef spade::ctrs::array<real_t, 3> v3d;
    
    spade::parallel::mpi_t group(&argc, &argv);
    
    spade::viscous_laws::constant_viscosity_t<real_t> visc_law(1.0);
    visc_law.prandtl = 1.0;
    
    spade::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.14285714285705;
    air.gamma = 1.4;
    
    spade::navier_stokes_mms::cns_perfectgas_mms_t mms(air, visc_law);
    auto mms_test_func = [&](const v3d& xyz) -> spade::fluid_state::prim_t<real_t>
    {
        return mms.test_fcn(xyz);
    };
    
    auto mms_visc_func = [&](const v3d& xyz) -> spade::ctrs::array<real_t,5>
    {
        return mms.visc_rhs(xyz);
    };
    
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) = -0.5;
    bounds.max(0) =  0.5;
    bounds.min(1) =  0.5;
    bounds.max(1) =  1.5;
    bounds.min(2) = -0.5;
    bounds.max(2) =  0.5;
    
    spade::coords::integrated_tanh_1D<real_t> xc(bounds.min(0), bounds.max(0), 0.1, 1.3);
    spade::coords::integrated_tanh_1D<real_t> yc(bounds.min(1), bounds.max(1), 0.1, 1.3);
    spade::coords::integrated_tanh_1D<real_t> zc(bounds.min(2), bounds.max(2), 0.1, 1.3);
    // spade::coords::quad_1D<real_t> xc;
    // spade::coords::quad_1D<real_t> yc;
    // spade::coords::quad_1D<real_t> zc;
    // spade::coords::scaled_coord_1D<real_t> xc(2.0);
    // spade::coords::scaled_coord_1D<real_t> yc(2.0);
    // spade::coords::scaled_coord_1D<real_t> zc(2.0);
    // spade::coords::identity_1D<real_t> xc;
    // spade::coords::identity_1D<real_t> yc;
    // spade::coords::identity_1D<real_t> zc;
    spade::coords::diagonal_coords coords(xc, yc, zc);
    // spade::coords::cyl_coords<real_t> coords;
    
    // spade::coords::identity<real_t> coords;
    spade::ctrs::array<std::size_t, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<std::size_t, 3> exchange_cells(2, 2, 2);
    
    spade::viscous::visc_lr visc_scheme(visc_law);
    
    std::vector<int> numcells = {8, 12, 16, 24, 30, 34, 38};
    
    std::vector<real_t> err_0_linf(numcells.size());
    std::vector<real_t> err_1_linf(numcells.size());
    std::vector<real_t> err_2_linf(numcells.size());
    std::vector<real_t> err_3_linf(numcells.size());
    std::vector<real_t> err_4_linf(numcells.size());
    
    std::vector<real_t> err_0_l2  (numcells.size());
    std::vector<real_t> err_1_l2  (numcells.size());
    std::vector<real_t> err_2_l2  (numcells.size());
    std::vector<real_t> err_3_l2  (numcells.size());
    std::vector<real_t> err_4_l2  (numcells.size());
    std::vector<real_t> dxs;    
    
    for (int n = 0; n < numcells.size(); ++n)
    {
        spade::ctrs::array<std::size_t, 3> cells_in_block(numcells[n], numcells[n], numcells[n]);
    
        spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
        dxs.push_back(spade::utils::min(grid.get_dx(0), grid.get_dx(1), grid.get_dx(2)));
    
        std::size_t num_dof = grid.get_num_interior_cells();
        num_dof = group.sum(num_dof);
    
        if (group.isroot()) print("Grid level / cells:           ", n, "/", num_dof);
    
        prim_t fill1 = 0.0;
        flux_t fill2 = 0.0;
        spade::grid::grid_array prim    (grid, fill1);
        spade::grid::grid_array rhs     (grid, fill2);
        spade::grid::grid_array rhs_test(grid, fill2);
    
        spade::algs::fill_array(prim,     mms_test_func, spade::grid::include_exchanges);
        spade::algs::fill_array(rhs_test, mms_visc_func, spade::grid::include_exchanges);        
    
        spade::pde_algs::flux_div(prim, rhs, visc_scheme);
    
        bool output = true;
        if (output)
        {
            std::string filename = "rhs" + std::to_string(n);
            std::string out_file = spade::io::output_vtk("output", filename, grid, rhs);
            filename = "ana" + std::to_string(n);
            out_file = spade::io::output_vtk("output", filename, grid, rhs_test);
        }
        rhs -= rhs_test;
    
        spade::reduce_ops::reduce_sum<real_t> sum;
        spade::reduce_ops::reduce_max<real_t> max;
    
        err_1_l2[n]   = sqrt(spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[1]*rhs_arr[1];}, sum)/num_dof);
        err_2_l2[n]   = sqrt(spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[2]*rhs_arr[2];}, sum)/num_dof);
        err_3_l2[n]   = sqrt(spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[3]*rhs_arr[3];}, sum)/num_dof);
        err_4_l2[n]   = sqrt(spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[4]*rhs_arr[4];}, sum)/num_dof);
    
        err_1_linf[n] = spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[1]);}, max);
        err_2_linf[n] = spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[2]);}, max);
        err_3_linf[n] = spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[3]);}, max);
        err_4_linf[n] = spade::algs::transform_reduce(rhs, [&](const spade::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[4]);}, max);
    
        if (group.isroot()) print("L2:  ", err_1_l2  [n], err_2_l2  [n], err_3_l2  [n], err_4_l2  [n]);
        if (group.isroot()) print("LInf:", err_1_linf[n], err_2_linf[n], err_3_linf[n], err_4_linf[n]);
        if (group.isroot()) print();
    }
    auto err_report = [](const std::string& title, const std::vector<real_t>& errs, const std::vector<real_t>& dx) -> void
    {
        std::vector<real_t> orders(errs.size()-1);
        for (int i = 0; i < orders.size(); ++i)
        {
            orders[i] = (log(errs[i+1])-log(errs[i]))/(log(dx[i+1])-log(dx[i]));
        }
        print(">>", title);
        for (auto r: orders)
        {
            print(r);
        }
    };
    if (group.isroot())
    {
        print("Error report:");
        print("=========================================");
        err_report("Energy     (conv, L2)", err_1_l2, dxs);
        err_report("X-Momentum (conv, L2)", err_2_l2, dxs);
        err_report("Y-Momentum (conv, L2)", err_3_l2, dxs);
        err_report("Z-Momentum (conv, L2)", err_4_l2, dxs);
        print("=========================================");
        err_report("Energy     (conv, Linf)", err_1_linf, dxs);
        err_report("X-Momentum (conv, Linf)", err_2_linf, dxs);
        err_report("Y-Momentum (conv, Linf)", err_3_linf, dxs);
        err_report("Z-Momentum (conv, Linf)", err_4_linf, dxs);
        print("=========================================");
        print("Happy computing.");
    }
    
    return 0;
}