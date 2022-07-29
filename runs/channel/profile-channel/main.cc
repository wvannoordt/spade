#include "PTL.h"
#include "spade.h"
#include "local_types.h"
#include "dns_filter.h"
#include "prof_t.h"
#include "calc_u_bulk.h"

int main(int argc, char** argv)
{
    
    for (auto pp: range(6,5)*range(0,5))
    {
        print(pp[0], pp[1]);
    }
    return 0;
/*    
    spade::parallel::mpi_t group(&argc, &argv);
    
    std::string input_filename = "input.ptl";
    PTL::PropertyTree input;
    // input.Read(input_filename);
    // int filt_x  = input["filt_x"];
    // int filt_y  = input["filt_y"];
    // int filt_z  = input["filt_z"];
    // bool output = input["output"];
    int filt_x  = 13;
    int filt_y  = 13;
    int filt_z  = 13;
    bool output = false;
    v3i filt(filt_x, filt_y, filt_z);
    
    spade::ctrs::array<int, 3> num_blocks(8, 8, 8);
    spade::ctrs::array<int, 3> cells_in_block(48, 48, 48);
    spade::ctrs::array<int, 3> cells_in_block_coarse;
    for (auto i: range(0,3)) cells_in_block_coarse[i] = cells_in_block[i]/filt[i];
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::ctrs::array<int, 3> exchange_cells_filt(filt[0]/2, filt[1]/2, filt[2]/2);
    spade::bound_box_t<real_t, 3> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*spade::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*spade::consts::pi*delta;
    prim_t fill1 = 0;
    spade::coords::identity<real_t> coords;
    
    spade::grid::cartesian_grid_t grid     (num_blocks, cells_in_block,        exchange_cells,      bounds, coords, group);
    spade::grid::cartesian_grid_t grid_filt(num_blocks, cells_in_block,        exchange_cells_filt, bounds, coords, group);
    // spade::grid::cartesian_grid_t grid_crse(num_blocks, cells_in_block_coarse, exchange_cells,      bounds, coords, group);
    
    spade::grid::grid_array prim(grid, fill1);
    spade::grid::grid_array prim_o(grid_filt, fill1);
    spade::grid::grid_array prim_i(grid_filt, fill1);
    
    
    
    spade::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    spade::fluid_state::perfect_gas_t<real_t> air;
    air.R = 287.15;
    air.gamma = 1.4;
    
    const real_t p0         = 30.0;
    const real_t t0         = 0.1;
    const real_t u0         = 69.54;
    const real_t mu         = visc_law.get_visc();
    const real_t rho        = p0/(air.R*t0);
    const real_t u_tau      = re_tau*mu/(rho*delta);
    const real_t force_term = rho*u_tau*u_tau/delta;
    const real_t du         = 3.0;
    
    
    const int nyi = grid.get_num_cells(1)*grid.get_num_blocks(1);
    const int nyo = nyi/filt[1];
    const real_t dyi = 2*delta/nyi;
    const real_t dyo = 2*delta/nyo;
    const int ny = nyi;
    
    std::vector<profr_t*> reg;
    
    profr_t y_f    (ny+1, 0.0, "y_f",     reg); // 1
    profr_t ui_f   (ny+1, 0.0, "ui_f",    reg); // 2
    profr_t uo_f   (ny+1, 0.0, "uo_f",    reg); // 3
    profr_t vi_f   (ny+1, 0.0, "vi_f",    reg); // 4
    profr_t vo_f   (ny+1, 0.0, "vo_f",    reg); // 5
    profr_t wi_f   (ny+1, 0.0, "wi_f",    reg); // 6
    profr_t wo_f   (ny+1, 0.0, "wo_f",    reg); // 7
    profr_t ui2_f  (ny+1, 0.0, "ui2_f",   reg); // 8
    profr_t uo2_f  (ny+1, 0.0, "uo2_f",   reg); // 9
    profr_t vi2_f  (ny+1, 0.0, "vi2_f",   reg); // 10
    profr_t vo2_f  (ny+1, 0.0, "vo2_f",   reg); // 11
    profr_t wi2_f  (ny+1, 0.0, "wi2_f",   reg); // 12
    profr_t wo2_f  (ny+1, 0.0, "wo2_f",   reg); // 13
    profr_t uivi_f (ny+1, 0.0, "uivi_f",  reg); // 14
    profr_t uovo_f (ny+1, 0.0, "uovo_f",  reg); // 15
    profr_t uiuo_f (ny+1, 0.0, "uiuo_f",  reg); // 16
    profr_t vivo_f (ny+1, 0.0, "vivo_f",  reg); // 17
    profr_t wiwo_f (ny+1, 0.0, "wiwo_f",  reg); // 18
    profr_t uivo_f (ny+1, 0.0, "uivo_f",  reg); // 19
    profr_t viuo_f (ny+1, 0.0, "viuo_f",  reg); // 20
    profr_t duidy_f(ny+1, 0.0, "duidy_f", reg); // 21
    profr_t duody_f(ny+1, 0.0, "duody_f", reg); // 22
    
    profr_t y    (ny+1, 0.0, "y",    reg); // 23
    profr_t ui   (ny+1, 0.0, "ui",   reg); // 24
    profr_t uo   (ny+1, 0.0, "uo",   reg); // 25
    profr_t vi   (ny+1, 0.0, "vi",   reg); // 26
    profr_t vo   (ny+1, 0.0, "vo",   reg); // 27
    profr_t wi   (ny+1, 0.0, "wi",   reg); // 28
    profr_t wo   (ny+1, 0.0, "wo",   reg); // 29
    profr_t ui2  (ny+1, 0.0, "ui2",  reg); // 30
    profr_t uo2  (ny+1, 0.0, "uo2",  reg); // 31
    profr_t vi2  (ny+1, 0.0, "vi2",  reg); // 32
    profr_t vo2  (ny+1, 0.0, "vo2",  reg); // 33
    profr_t wi2  (ny+1, 0.0, "wi2",  reg); // 34
    profr_t wo2  (ny+1, 0.0, "wo2",  reg); // 35
    profr_t uivi (ny+1, 0.0, "uivi", reg); // 36
    profr_t uovo (ny+1, 0.0, "uovo", reg); // 37
    profr_t uiuo (ny+1, 0.0, "uiuo", reg); // 38
    profr_t vivo (ny+1, 0.0, "vivo", reg); // 39
    profr_t wiwo (ny+1, 0.0, "wiwo", reg); // 40
    profr_t uivo (ny+1, 0.0, "uivo", reg); // 41
    profr_t viuo (ny+1, 0.0, "viuo", reg); // 42
    profr_t p1   (ny+1, 0.0, "p1",   reg); // 43
    profr_t p2   (ny+1, 0.0, "p2",   reg); // 44
    
    std::vector<std::string> names;
    for (int i = 1; i < argc; i++) names.push_back(std::string(argv[i]));
    std::ofstream ub_out("ub.dat");
    int ct = 0;
    for (auto& p: names)
    {
        for (auto symmetry_index: range(0,2))
        {
            if (group.isroot())
            {
                print("File:", p, "  symmetry:", symmetry_index);
            }
            if (!std::filesystem::exists(p))
            {
                if (group.isroot()) print("The following file does not exsist:", p);
                group.sync();
                return 15;
            }
            spade::io::binary_read(p, prim);
            m3r symmetry_jacobian;
            switch (symmetry_index)
            {
                case 0: { symmetry_jacobian = m3r({{ 1.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0 },{ 0.0, 0.0, 1.0 }}); break; }
                case 1: { symmetry_jacobian = m3r({{ 1.0, 0.0, 0.0 },{ 0.0,-1.0, 0.0 },{ 0.0, 0.0, 1.0 }}); break; }
                case 2: { symmetry_jacobian = m3r({{ 1.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0 },{ 0.0, 0.0,-1.0 }}); break; }
                case 3: { symmetry_jacobian = m3r({{ 1.0, 0.0, 0.0 },{ 0.0,-1.0, 0.0 },{ 0.0, 0.0,-1.0 }}); break; }
            }
            spade::algs::transform_inplace(prim, [&](const prim_t& q) -> prim_t
            {
                prim_t q_out;
                q_out.p() = q.p();
                q_out.T() = q.T();
                v3r u_vec(q.u(), q.v(), q.w());
                u_vec = symmetry_jacobian * u_vec;
                q_out.u() = u_vec[0];
                q_out.v() = u_vec[1];
                q_out.w() = u_vec[2];
                return q_out;
            });
            v3r e_y(0,1,0);
            v3r e_y_sym = symmetry_jacobian*e_y;
            postprocessing::copy_field(prim, prim_i);
            grid_filt.exchange_array(prim_i);
    
            const real_t ub = calc_u_bulk(prim, air);
            if (group.isroot()) ub_out << ub << std::endl;
            postprocessing::noslip(filt[1]/2, prim_i);
            postprocessing::noslip(filt[1]/2, prim_o);
            postprocessing::spatial_filter(filt, prim_i, prim_o);
            prim_i -= prim_o;
            grid_filt.exchange_array(prim_i);
            grid_filt.exchange_array(prim_o);
            if (output)
            {
                postprocessing::copy_field(prim_i, prim);
                if (group.isroot()) print("Outputting...");
                spade::io::output_vtk("output", "q_i", prim);
                postprocessing::copy_field(prim_o, prim);
                spade::io::output_vtk("output", "q_o", prim);
                if (group.isroot()) print("Done. Exiting.");
                group.sync();
                return 0;
            }
    
            postprocessing::extract_profile(symmetry_jacobian, y,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return x[1];});
            postprocessing::extract_profile(symmetry_jacobian, ui,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u();});
            postprocessing::extract_profile(symmetry_jacobian, uo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u();});
            postprocessing::extract_profile(symmetry_jacobian, vi,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v();});
            postprocessing::extract_profile(symmetry_jacobian, vo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.v();});
            postprocessing::extract_profile(symmetry_jacobian, wi,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w();});
            postprocessing::extract_profile(symmetry_jacobian, wo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.w();});
            postprocessing::extract_profile(symmetry_jacobian, ui2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_i.u();});
            postprocessing::extract_profile(symmetry_jacobian, uo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_o.u();});
            postprocessing::extract_profile(symmetry_jacobian, vi2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v()*q_i.v();});
            postprocessing::extract_profile(symmetry_jacobian, vo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.v()*q_o.v();});
            postprocessing::extract_profile(symmetry_jacobian, wi2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w()*q_i.w();});
            postprocessing::extract_profile(symmetry_jacobian, wo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.w()*q_o.w();});
            postprocessing::extract_profile(symmetry_jacobian, uivi, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_i.v();});
            postprocessing::extract_profile(symmetry_jacobian, uovo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_o.v();});
            postprocessing::extract_profile(symmetry_jacobian, uiuo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_o.u();});
            postprocessing::extract_profile(symmetry_jacobian, vivo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v()*q_o.v();});
            postprocessing::extract_profile(symmetry_jacobian, wiwo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w()*q_o.w();});
            postprocessing::extract_profile(symmetry_jacobian, uivo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_o.v();});
            postprocessing::extract_profile(symmetry_jacobian, viuo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_i.v();});
            postprocessing::extract_profile(symmetry_jacobian, p1,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return (q_i.p()+q_o.p());});
            postprocessing::extract_profile(symmetry_jacobian, p2,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return (q_i.p()+q_o.p())*(q_i.p()+q_o.p());});
    
            auto sqr = [](const real_t& x) -> real_t {return x*x;};
            postprocessing::extract_profile_f(symmetry_jacobian, y_f,     prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return x[1]; });
            postprocessing::extract_profile_f(symmetry_jacobian, ui_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u()); });
            postprocessing::extract_profile_f(symmetry_jacobian, uo_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.u()+q_o_R.u()); });
            postprocessing::extract_profile_f(symmetry_jacobian, vi_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.v()+q_i_R.v()); });
            postprocessing::extract_profile_f(symmetry_jacobian, vo_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.v()+q_o_R.v()); });
            postprocessing::extract_profile_f(symmetry_jacobian, wi_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.w()+q_i_R.w()); });
            postprocessing::extract_profile_f(symmetry_jacobian, wo_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.w()+q_o_R.w()); });
            postprocessing::extract_profile_f(symmetry_jacobian, ui2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_i_L.u()+q_i_R.u())); });
            postprocessing::extract_profile_f(symmetry_jacobian, uo2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_o_L.u()+q_o_R.u())); });
            postprocessing::extract_profile_f(symmetry_jacobian, vi2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_i_L.v()+q_i_R.v())); });
            postprocessing::extract_profile_f(symmetry_jacobian, vo2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_o_L.v()+q_o_R.v())); });
            postprocessing::extract_profile_f(symmetry_jacobian, wi2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_i_L.w()+q_i_R.w())); });
            postprocessing::extract_profile_f(symmetry_jacobian, wo2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_o_L.w()+q_o_R.w())); });
            postprocessing::extract_profile_f(symmetry_jacobian, uivi_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u())*0.5*(q_i_L.v()+q_i_R.v()); });
            postprocessing::extract_profile_f(symmetry_jacobian, uovo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.u()+q_o_R.u())*0.5*(q_o_L.v()+q_o_R.v()); });
            postprocessing::extract_profile_f(symmetry_jacobian, uiuo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u())*0.5*(q_o_L.u()+q_o_R.u()); });
            postprocessing::extract_profile_f(symmetry_jacobian, vivo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.v()+q_i_R.v())*0.5*(q_o_L.v()+q_o_R.v()); });
            postprocessing::extract_profile_f(symmetry_jacobian, wiwo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.w()+q_i_R.w())*0.5*(q_o_L.w()+q_o_R.w()); });
            postprocessing::extract_profile_f(symmetry_jacobian, uivo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u())*0.5*(q_o_L.v()+q_o_R.v()); });
            postprocessing::extract_profile_f(symmetry_jacobian, viuo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.v()+q_i_R.v())*0.5*(q_o_L.u()+q_o_R.u()); });
            postprocessing::extract_profile_f(symmetry_jacobian, duody_f, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return e_y_sym[1]*(q_o_R.u()-q_o_L.u())/dyo;});
            postprocessing::extract_profile_f(symmetry_jacobian, duidy_f, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return e_y_sym[1]*(q_i_R.u()-q_i_L.u())/dyi;});
    
            for (auto p:reg) p->aggregate();
        }
    }
    bool output_names = false;
    if (group.isroot())
    {
        std::filesystem::create_directory("profiles");
        std::ofstream myfile("profiles/pfs.dat");
        if (output_names)
        {
            for (int n = 0; n < reg.size(); ++n) myfile << reg[n]->name << ((n<(reg.size()-1))?",":"");
            myfile << "\n";
        }
        for (int k = 0; k < reg[0]->avg.size(); ++k)
        {
            for (int n = 0; n < reg.size(); ++n)
            {
                auto& vec = reg[n]->avg;
                myfile << vec[k] << ((n<(reg.size()-1))?",":"");
            }
            myfile << "\n";
        }
    }
    if (group.isroot()) print("done!");
    return 0;*/
}
