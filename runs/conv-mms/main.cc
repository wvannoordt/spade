#include "cvdf.h"

int main(int argc, char** argv)
{
    typedef double real_t;
    typedef cvdf::ctrs::array<real_t, 3> v3d;
    
    
    
    
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.14285714285705;
    air.gamma = 1.4;
    
    cvdf::navier_stokes_mms::cns_perfectgas_mms_t mms(air, visc_law);
    auto mms_test_func = [&](const v3d& xyz) -> cvdf::fluid_state::prim_t<real_t>
    {
        return mms.test_fcn(xyz);
    };
    
    auto mms_conv_func = [&](const v3d& xyz) -> cvdf::ctrs::array<real_t,5>
    {
        return mms.conv_rhs(xyz);
    };
    
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    bounds.min(0) = -0.5;
    bounds.max(0) =  0.5;
    bounds.min(1) =  1.5;
    bounds.max(1) =  2.5;
    bounds.min(2) = -0.5;
    bounds.max(2) =  0.5;
    
    // cvdf::coords::integrated_tanh_1D<real_t> xc(bounds.min(0), bounds.max(0), 0.1, 1.3);
    cvdf::coords::integrated_tanh_1D<real_t> yc(bounds.min(1), bounds.max(1), 0.1, 1.3);
    // cvdf::coords::integrated_tanh_1D<real_t> zc(bounds.min(2), bounds.max(2), 0.1, 1.3);
    
    cvdf::coords::identity_1D<real_t> xc;
    // cvdf::coords::identity_1D<real_t> yc;
    cvdf::coords::identity_1D<real_t> zc;
    cvdf::coords::diagonal_coords coords(xc, yc, zc);
    // cvdf::coords::cyl_coords<real_t> coords;
    
    // cvdf::coords::identity<real_t> coords;
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> num_blocks(2, 2, 2);
    cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    
    cvdf::convective::totani_lr tscheme(air);
    
    // std::vector<int> numcells = {24, 26};
    std::vector<int> numcells = {8, 12, 16, 24};
    
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
        cvdf::ctrs::array<std::size_t, cvdf::cvdf_dim> cells_in_block(numcells[n], numcells[n], numcells[n]);
        
        cvdf::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
        
        dxs.push_back(cvdf::utils::min(grid.get_dx(0), grid.get_dx(1), grid.get_dx(2)));
        
        std::size_t num_dof = grid.get_num_interior_cells();
        num_dof = group.sum(num_dof);
        
        if (group.isroot()) print("Grid level / cells:           ", n, "/", num_dof);
        
        cvdf::grid::grid_array prim    (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
        cvdf::grid::grid_array rhs     (grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
        cvdf::grid::grid_array rhs_test(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
        
        cvdf::algs::fill_array(prim,     mms_test_func, cvdf::grid::include_exchanges);
        cvdf::algs::fill_array(rhs_test, mms_conv_func, cvdf::grid::include_exchanges);        
        
        cvdf::flux_algs::flux_lr_diff(prim, rhs, tscheme);
        
        bool output = true;
        if (output)
        {
            std::string filename = "rhs" + std::to_string(n);
            std::string out_file = cvdf::output::output_vtk("output", filename, grid, rhs);
        }
        rhs -= rhs_test;
        
        cvdf::reduce_ops::reduce_sum<real_t> sum;
        cvdf::reduce_ops::reduce_max<real_t> max;
        
        err_0_l2[n]   = sqrt(cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[0]*rhs_arr[0];}, sum)/num_dof);
        err_1_l2[n]   = sqrt(cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[1]*rhs_arr[1];}, sum)/num_dof);
        err_2_l2[n]   = sqrt(cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[2]*rhs_arr[2];}, sum)/num_dof);
        err_3_l2[n]   = sqrt(cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[3]*rhs_arr[3];}, sum)/num_dof);
        err_4_l2[n]   = sqrt(cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return rhs_arr[4]*rhs_arr[4];}, sum)/num_dof);
        
        err_0_linf[n] = cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[0]);}, max);
        err_1_linf[n] = cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[1]);}, max);
        err_2_linf[n] = cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[2]);}, max);
        err_3_linf[n] = cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[3]);}, max);
        err_4_linf[n] = cvdf::algs::transform_reduce(rhs, [&](const cvdf::ctrs::array<real_t, 5>& rhs_arr) -> real_t {return std::abs(rhs_arr[4]);}, max);
        
        if (group.isroot()) print("L2:  ", err_0_l2  [n], err_1_l2  [n], err_2_l2  [n], err_3_l2  [n], err_4_l2  [n]);
        if (group.isroot()) print("LInf:", err_0_linf[n], err_1_linf[n], err_2_linf[n], err_3_linf[n], err_4_linf[n]);
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
        err_report("Continuity (conv, L2)", err_0_l2, dxs);
        err_report("Energy     (conv, L2)", err_1_l2, dxs);
        err_report("X-Momentum (conv, L2)", err_2_l2, dxs);
        err_report("Z-Momentum (conv, L2)", err_3_l2, dxs);
        err_report("Y-Momentum (conv, L2)", err_4_l2, dxs);
        print("=========================================");
        err_report("Continuity (conv, Linf)", err_0_linf, dxs);
        err_report("Energy     (conv, Linf)", err_1_linf, dxs);
        err_report("X-Momentum (conv, Linf)", err_2_linf, dxs);
        err_report("Z-Momentum (conv, Linf)", err_3_linf, dxs);
        err_report("Y-Momentum (conv, Linf)", err_4_linf, dxs);
        print("=========================================");
        print("Happy computing.");
    }
    
    return 0;
}