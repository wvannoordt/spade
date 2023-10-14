#include <stdio.h>
#include <vector>

#include "spade.h"

using real_t = double;
using flux_t = spade::fluid_state::flux_t<real_t>;
using prim_t = spade::fluid_state::prim_t<real_t>;
using cons_t = spade::fluid_state::cons_t<real_t>;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
        
    spade::ctrs::array<int, 3>    num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3>    cells_in_block(16, 16, 16);
    spade::ctrs::array<int, 3>    exchange_cells(2, 2, 2);
    // spade::bound_box_t<real_t, 3> bounds(-0.2, 1.2, -0.4, 0.4, -0.4, 0.4);
    // spade::bound_box_t<real_t, 3> bounds(-0.5, 1.5, -0.5, 0.5, 0.0, 1.0);
    spade::bound_box_t<real_t, 3> bounds(-2.2, 2.2, -2.2, 2.2, -2.2, 2.2);
    
    spade::coords::identity<real_t> coords;
    
    spade::ctrs::array<bool, 3> periodic = false;
    spade::amr::amr_blocks_t blocks(num_blocks, bounds);
    spade::grid::cartesian_grid_t grid(cells_in_block, exchange_cells, blocks, coords, group);
    
    spade::geom::vtk_geom_t<3> geom;
    spade::geom::read_vtk_geom("sphere.vtk", geom);
    
    using bvh_t = spade::geom::bvh_t<2, real_t>;
    bvh_t b0;
    
    spade::bound_box_t<real_t, 2> bnd;
    bnd.min(0) = geom.bbox_inflated.min(0);
    bnd.max(0) = geom.bbox_inflated.max(0);
    bnd.min(1) = geom.bbox_inflated.min(1);
    bnd.max(1) = geom.bbox_inflated.max(1);
    
    const auto check = [&](const std::size_t& i, const auto& ibnd)
    {
        spade::bound_box_t<real_t, 3> ext_bx;
        ext_bx.min(0) = ibnd.min(0);
        ext_bx.max(0) = ibnd.max(0);
        ext_bx.min(1) = ibnd.min(1);
        ext_bx.max(1) = ibnd.max(1);
        ext_bx.min(2) = 0.0;
        ext_bx.max(2) = 1.0;
        
        auto ps = geom.faces[i];
        auto p0 = geom.points[ps[0]];
        auto p1 = geom.points[ps[1]];
        auto p2 = geom.points[ps[2]];
        
        p0[2] = -0.01;
        p1[2] = 1.01;
        p2[2] = 1.01;
        
        return spade::geom::primitives::box_tri_intersect(p0, p1, p2, ext_bx);
    };
    
    bool do_refine = true;
    if (do_refine)
    {
        spade::timing::scoped_tmr_t t0("refine");
        int maxlevel = 2;
        {
            while (true)
            {
                const auto bndy_intersect = [&](const auto& lb)
                {
                    if (grid.get_blocks().get_amr_node(lb).level[0] >= maxlevel) return false;
                    const auto bnd = grid.get_bounding_box(lb);
                    return geom.box_contains_boundary(bnd);
                };
                
                auto rblks = grid.select_blocks(bndy_intersect, spade::partition::global);
                if (rblks.size() == 0) break;
                grid.refine_blocks(rblks);
            }
        }
    }
    if (group.isroot())
    {
        std::size_t npt = 1;
        npt *= cells_in_block[0];
        npt *= cells_in_block[1];
        npt *= cells_in_block[2];
        npt *= grid.get_num_global_blocks();
        print("num points:", npt);
        print("num tri:   ", geom.faces.size());
    }
    
    const auto ghosts = spade::ibm::compute_ghosts(grid, geom);
    spade::io::output_vtk("pts.vtk", ghosts.boundary_points.data(spade::device::cpu));
    spade::io::output_vtk("cls.vtk", ghosts.closest_points.data(spade::device::cpu));
    
    std::vector<spade::coords::point_t<real_t>> ips;
    int ng = ghosts.indices.data(spade::device::cpu).size();
    for (int i = 0; i < ng; ++i)
    {
        const auto& icell = ghosts.indices[i];
        const auto xg = grid.get_comp_coords(icell);
        const auto bp = ghosts.closest_points[i];
        auto nv = bp - xg;
        nv = nv / spade::ctrs::array_norm(nv);
        auto ipt = xg + 0.05*nv;
        ips.push_back(ipt);
    }
    
    spade::io::output_vtk("ips.vtk", ips);
    
    using prim_t = spade::fluid_state::prim_t<real_t>;
    prim_t finit = 0.0;
    spade::grid::grid_array phi(grid, finit, spade::device::best);
    auto phi_img = phi.image();
    int ct = 0;
    for (const auto& ig: ghosts.indices.data(spade::device::cpu))
    {
        auto ff         = phi_img.get_elem(ig);
        const auto idir = ghosts.directions[ct];
        
        if (idir == 0) ff.u() += 1.0;
        if (idir == 1) ff.v() += 1.0;
        if (idir == 2) ff.w() += 1.0;
        
        phi_img.set_elem(ig, ff);
        ++ct;
    }
    spade::io::output_vtk("output", "phi", phi);
    return 0;
}
