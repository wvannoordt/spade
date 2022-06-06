#include "cvdf.h"
#include "local_types.h"
#include "dns_filter.h"
#include "prof_t.h"

int main(int argc, char** argv)
{
    cvdf::parallel::mpi_t group(&argc, &argv);
    
    
    cvdf::ctrs::array<int, cvdf::cvdf_dim> num_blocks(8, 8, 8);
    cvdf::ctrs::array<int, cvdf::cvdf_dim> cells_in_block(48, 48, 48);
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
    
    cvdf::coords::identity<real_t> coords;
    
    cvdf::grid::cartesian_grid_t grid     (num_blocks, cells_in_block, exchange_cells,      bounds, coords, group);
    cvdf::grid::cartesian_grid_t grid_filt(num_blocks, cells_in_block, exchange_cells_filt, bounds, coords, group);
    
    cvdf::grid::grid_array prim(grid, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array prim_o(grid_filt, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    cvdf::grid::grid_array prim_i(grid_filt, 0.0, cvdf::dims::static_dims<5>(), cvdf::dims::singleton_dim());
    
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
    
    const int ny = grid.get_num_cells(1)*grid.get_num_blocks(1);
    
    std::vector<profr_t*> reg;
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
    profr_t uiuo (ny, 0.0, "uiuo", reg);
    profr_t vivo (ny, 0.0, "vivo", reg);
    profr_t wiwo (ny, 0.0, "wiwo", reg);
    profr_t uivo (ny, 0.0, "uivo", reg);
    profr_t viuo (ny, 0.0, "viuo", reg);
    
    std::vector<std::string> names;
    for (int i = 1; i < argc; i++) names.push_back(std::string(argv[i]));
    
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
        postprocessing::dns_filter(prim_i, prim_o);
        prim_i -= prim_o;
        postprocessing::extract_profile(y,    prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return x[1];});
        postprocessing::extract_profile(ui,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u();});
        postprocessing::extract_profile(uo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u();});
        postprocessing::extract_profile(vi,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v();});
        postprocessing::extract_profile(vo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.v();});
        postprocessing::extract_profile(wi,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u();});
        postprocessing::extract_profile(wo,   prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u();});
        postprocessing::extract_profile(ui2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_i.u();});
        postprocessing::extract_profile(uo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_o.u();});
        postprocessing::extract_profile(vi2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v()*q_i.v();});
        postprocessing::extract_profile(vo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.v()*q_o.v();});
        postprocessing::extract_profile(wi2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w()*q_i.w();});
        postprocessing::extract_profile(wo2,  prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.w()*q_o.w();});
        postprocessing::extract_profile(uiuo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_o.u();});
        postprocessing::extract_profile(vivo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.v()*q_o.v();});
        postprocessing::extract_profile(wiwo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.w()*q_o.w();});
        postprocessing::extract_profile(uivo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_i.u()*q_o.v();});
        postprocessing::extract_profile(viuo, prim_o, prim_i, [&](const v3d& x, const prim_t& q_o, const prim_t& q_i) -> real_t {return q_o.u()*q_i.v();});
        for (auto p:reg) p->aggregate();
    }
    if (group.isroot())
    {
        std::ofstream myfile("profs.dat");
        for (int n = 0; n < reg.size(); ++n) myfile << reg[n]->name << ((n<(reg.size()-1))?",":"");
        myfile << "\n";
        for (int k = 0; k < reg[0]->avg.size(); ++k)
        {
            for (int n = 0; n < reg.size(); ++n)
            {
                auto& vec = reg[k]->avg;
                myfile << vec[k] << ((n<(reg.size()-1))?",":"");
            }
            myfile << "\n";
        }
    }
    return 0;
}
