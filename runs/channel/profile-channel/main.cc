#include "cvdf.h"
#include "local_types.h"
#include "dns_filter.h"
#include "prof_t.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    v3i filt(8, 16, 6);
    //v3i filt(1, 1, 1);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> num_blocks(8, 8, 8);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block(48, 48, 48);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block_coarse;
    for (auto i: range(0,3)) cells_in_block_coarse[i] = cells_in_block[i]/filt[i];
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> exchange_cells_filt(8, 8, 8);
    cvdf::bound_box_t<real_t, cvdf::cvdf_dim> bounds;
    const real_t re_tau = 180.0;
    const real_t delta = 1.0;
    bounds.min(0) =  0.0;
    bounds.max(0) =  4.0*cvdf::consts::pi*delta;
    bounds.min(1) = -delta;
    bounds.max(1) =  delta;
    bounds.min(2) =  0.0;
    bounds.max(2) =  2*cvdf::consts::pi*delta;
    prim_t fill1 = 0;
    cvdf::coords::identity<real_t> coords;
    
    cvdf::grid::cartesian_grid_t grid     (num_blocks, cells_in_block,        exchange_cells,      bounds, coords, group);
    cvdf::grid::cartesian_grid_t grid_filt(num_blocks, cells_in_block,        exchange_cells_filt, bounds, coords, group);
    // cvdf::grid::cartesian_grid_t grid_crse(num_blocks, cells_in_block_coarse, exchange_cells,      bounds, coords, group);
    
    cvdf::grid::grid_array prim(grid, fill1);
    cvdf::grid::grid_array prim_o(grid_filt, fill1);
    cvdf::grid::grid_array prim_i(grid_filt, fill1);
    
    cvdf::viscous_laws::constant_viscosity_t<real_t> visc_law(1.85e-4);
    visc_law.prandtl = 0.72;
    
    cvdf::fluid_state::perfect_gas_t<real_t> air;
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
    
    profr_t y_f    (ny, 0.0, "y_f",     reg); // 1
    profr_t ui_f   (ny, 0.0, "ui_f",    reg); // 2
    profr_t uo_f   (ny, 0.0, "uo_f",    reg); // 3
    profr_t vi_f   (ny, 0.0, "vi_f",    reg); // 4
    profr_t vo_f   (ny, 0.0, "vo_f",    reg); // 5
    profr_t wi_f   (ny, 0.0, "wi_f",    reg); // 6
    profr_t wo_f   (ny, 0.0, "wo_f",    reg); // 7
    profr_t ui2_f  (ny, 0.0, "ui2_f",   reg); // 8
    profr_t uo2_f  (ny, 0.0, "uo2_f",   reg); // 9
    profr_t vi2_f  (ny, 0.0, "vi2_f",   reg); // 10
    profr_t vo2_f  (ny, 0.0, "vo2_f",   reg); // 11
    profr_t wi2_f  (ny, 0.0, "wi2_f",   reg); // 12
    profr_t wo2_f  (ny, 0.0, "wo2_f",   reg); // 13
    profr_t uivi_f (ny, 0.0, "uivi_f",  reg); // 14
    profr_t uovo_f (ny, 0.0, "uovo_f",  reg); // 15
    profr_t uiuo_f (ny, 0.0, "uiuo_f",  reg); // 16
    profr_t vivo_f (ny, 0.0, "vivo_f",  reg); // 17
    profr_t wiwo_f (ny, 0.0, "wiwo_f",  reg); // 18
    profr_t uivo_f (ny, 0.0, "uivo_f",  reg); // 19
    profr_t viuo_f (ny, 0.0, "viuo_f",  reg); // 20
    profr_t duidy_f(ny, 0.0, "duidy_f", reg); // 21
    profr_t duody_f(ny, 0.0, "duody_f", reg); // 22

    profr_t y    (ny, 0.0, "y",    reg);
    profr_t ui   (ny, 0.0, "ui",   reg);
    profr_t uo   (ny, 0.0, "uo",   reg);
    profr_t vi   (ny, 0.0, "vi",   reg);
    profr_t vo   (ny, 0.0, "vo",   reg);
    profr_t wi   (ny, 0.0, "wi",   reg);
    profr_t wo   (ny, 0.0, "wo",   reg);
    profr_t ui2  (ny, 0.0, "ui2",  reg);
    profr_t uo2  (ny, 0.0, "uo2",  reg);
    profr_t vi2  (ny, 0.0, "vi2",  reg);
    profr_t vo2  (ny, 0.0, "vo2",  reg);
    profr_t wi2  (ny, 0.0, "wi2",  reg);
    profr_t wo2  (ny, 0.0, "wo2",  reg);
    profr_t uivi (ny, 0.0, "uivi", reg);
    profr_t uovo (ny, 0.0, "uovo", reg);
    profr_t uiuo (ny, 0.0, "uiuo", reg);
    profr_t vivo (ny, 0.0, "vivo", reg);
    profr_t wiwo (ny, 0.0, "wiwo", reg);
    profr_t uivo (ny, 0.0, "uivo", reg);
    profr_t viuo (ny, 0.0, "viuo", reg);
    profr_t p1   (ny, 0.0, "p1",   reg);
    profr_t p2   (ny, 0.0, "p2",   reg);
    
    std::vector<std::string> names;
    for (int i = 1; i < argc; i++) names.push_back(std::string(argv[i]));
    
    bool output = false;
    int ct = 0;
    for (auto& p: names)
    {
        if (group.isroot()) print(p);
        if (!std::filesystem::exists(p))
        {
            if (group.isroot()) print("The following file does not exsist:", p);
            group.sync();
            return 15;
        }
        cvdf::io::binary_read(p, prim);
        postprocessing::copy_field(prim, prim_i);
        grid_filt.exchange_array(prim_i);
        postprocessing::dns_filter(filt, prim_i, prim_o);
	//        prim_i -= prim_o;
        grid_filt.exchange_array(prim_i);
        grid_filt.exchange_array(prim_o);
	postprocessing::noslip(filt[1]/2, prim_i);
        postprocessing::noslip(filt[1]/2, prim_o);
        if (output)
        {
            postprocessing::copy_field(prim_i, prim);
            cvdf::io::output_vtk("output", "q_i", prim);
            postprocessing::copy_field(prim_o, prim);
            cvdf::io::output_vtk("output", "q_o", prim);
        }
        postprocessing::extract_profile(y,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return x[1];});
        postprocessing::extract_profile(ui,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u();});
        postprocessing::extract_profile(uo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u();});
        postprocessing::extract_profile(vi,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v();});
        postprocessing::extract_profile(vo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.v();});
        postprocessing::extract_profile(wi,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w();});
        postprocessing::extract_profile(wo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.w();});
        postprocessing::extract_profile(ui2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_i.u();});
        postprocessing::extract_profile(uo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_o.u();});
        postprocessing::extract_profile(vi2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v()*q_i.v();});
        postprocessing::extract_profile(vo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.v()*q_o.v();});
        postprocessing::extract_profile(wi2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w()*q_i.w();});
        postprocessing::extract_profile(wo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.w()*q_o.w();});
        postprocessing::extract_profile(uivi, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_i.v();});
        postprocessing::extract_profile(uovo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_o.v();});
        postprocessing::extract_profile(uiuo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_o.u();});
        postprocessing::extract_profile(vivo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v()*q_o.v();});
        postprocessing::extract_profile(wiwo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w()*q_o.w();});
        postprocessing::extract_profile(uivo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_o.v();});
        postprocessing::extract_profile(viuo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_i.v();});
        postprocessing::extract_profile(p1,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return (q_i.p()+q_o.p());});
        postprocessing::extract_profile(p2,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return (q_i.p()+q_o.p())*(q_i.p()+q_o.p());});
        
        auto sqr = [](const real_t& x) -> real_t {return x*x;};
        postprocessing::extract_profile_f(y_f,     prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return x[1]; });
        postprocessing::extract_profile_f(ui_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u()); });
        postprocessing::extract_profile_f(uo_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.u()+q_o_R.u()); });
        postprocessing::extract_profile_f(vi_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.v()+q_i_R.v()); });
        postprocessing::extract_profile_f(vo_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.v()+q_o_R.v()); });
        postprocessing::extract_profile_f(wi_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.w()+q_i_R.w()); });
        postprocessing::extract_profile_f(wo_f,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.w()+q_o_R.w()); });
        postprocessing::extract_profile_f(ui2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_i_L.u()+q_i_R.u())); });
        postprocessing::extract_profile_f(uo2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_o_L.u()+q_o_R.u())); });
        postprocessing::extract_profile_f(vi2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_i_L.v()+q_i_R.v())); });
        postprocessing::extract_profile_f(vo2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_o_L.v()+q_o_R.v())); });
        postprocessing::extract_profile_f(wi2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_i_L.w()+q_i_R.w())); });
        postprocessing::extract_profile_f(wo2_f,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return sqr(0.5*(q_o_L.w()+q_o_R.w())); });
        postprocessing::extract_profile_f(uivi_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u())*0.5*(q_i_L.v()+q_i_R.v()); });
        postprocessing::extract_profile_f(uovo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_o_L.u()+q_o_R.u())*0.5*(q_o_L.v()+q_o_R.v()); });
        postprocessing::extract_profile_f(uiuo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u())*0.5*(q_o_L.u()+q_o_R.u()); });
        postprocessing::extract_profile_f(vivo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.v()+q_i_R.v())*0.5*(q_o_L.v()+q_o_R.v()); });
        postprocessing::extract_profile_f(wiwo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.w()+q_i_R.w())*0.5*(q_o_L.w()+q_o_R.w()); });
        postprocessing::extract_profile_f(uivo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.u()+q_i_R.u())*0.5*(q_o_L.v()+q_o_R.v()); });
        postprocessing::extract_profile_f(viuo_f,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return 0.5*(q_i_L.v()+q_i_R.v())*0.5*(q_o_L.u()+q_o_R.u()); });
        postprocessing::extract_profile_f(duody_f, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return (q_o_R.u()-q_o_L.u())/dyo;});
        postprocessing::extract_profile_f(duidy_f, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o_L, const prim_t& q_i_L, const prim_t& q_o_R, const prim_t& q_i_R) -> real_t { return (q_i_R.u()-q_i_L.u())/dyi;});
        
        for (auto p:reg) p->aggregate();
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
    return 0;
}
