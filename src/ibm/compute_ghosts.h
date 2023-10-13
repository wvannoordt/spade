#pragma once

#include "ibm/ghost_list.h"
#include "io/io_vtk.h"

namespace spade::ibm
{
    template <typename grid_t, typename geom_t>
    static ghost_list_t<typename grid_t::coord_type> compute_ghosts(const grid_t& grid, const geom_t& geom)
    {
        static_assert(std::same_as<typename grid_t::coord_sys_type, coords::identity<typename grid_t::coord_type>>, "ghosts currently only working for identity coordinates");
        using real_t = typename grid_t::coord_type;
        using pnt_t  = typename grid_t::coord_point_type;
        
        
        // Select the blocks that intersect the boundary
        const auto is_intersect = [&](const auto& lb_loc)
        {
            auto bnd_with_exch = grid.get_bounding_box(lb_loc);
            const auto dx = grid.get_dx(lb_loc);
            for (int d = 0; d < grid.dim(); ++d)
            {
                // bounding box includes exchange cells
                bnd_with_exch.min(d) -= grid.get_num_exchange(d)*dx[d];
                bnd_with_exch.max(d) += grid.get_num_exchange(d)*dx[d];
            }
            return geom.box_contains_boundary(bnd_with_exch);
        };
        const auto lbs = grid.select_blocks(is_intersect, partition::local);
        
        ghost_list_t<typename grid_t::coord_type> output;
        
        // We will only trace rays from the bounding box of the geometry
        const auto geom_bbx = geom.get_bounding_box().inflate(1.05);
        
        // Looping over boxes with boundary intersections
        for (const auto& lb: lbs)
        {
            // We will not consider xray intersection points outside of the domain
            const auto domain_boundary = grid.is_domain_boundary(lb);
            const auto dx = grid.get_dx(lb);
            for (int idir = 0; idir < grid.dim(); ++idir)
            {
                
                // Note that we include xray points inside the exchange cells in the tracing direction
                auto bnd_extd = grid.get_bounding_box(lb);
                if (!(domain_boundary.min(idir))) bnd_extd.min(idir) -= grid.get_num_exchange(idir)*dx[idir];
                if (!(domain_boundary.max(idir))) bnd_extd.max(idir) += grid.get_num_exchange(idir)*dx[idir];
                
                // Compute the indices that we will trace from
                bound_box_t<int, 3> ibound;
                for (int d = 0; d < grid.dim(); ++d)
                {
                    ibound.min(d) = 0;
                    ibound.max(d) = grid.get_num_cells(d);
                }
                ibound.min(idir) = 0;
                ibound.max(idir) = 1;
                
                using vec_t = ctrs::array<real_t, 3>;
                const auto loop_load = [&](const auto& ii)
                {
                    grid::cell_idx_t icell(ii[0], ii[1], ii[2], lb.value);
                    
                    // Note: this is performed in computational coordinates, so it is assumed that the geometry is in computational coordinates
                    auto x_comp  = grid.get_comp_coords(icell);
                    x_comp[idir] = geom_bbx.max(idir);
                    
                    // Now we trace the ray {x_comp, ray_vec}
                    
                    // We assume that we start out with a point outside
                    // the geometry (later modify this for internal flows)
                    bool from_exterior = false;
                    
                    // We will do this at every intersection point
                    int isect_count = 0;
                    const int sign = -1;
                    int sign_local = sign;
                    const auto on_intersection = [&](const auto& point)
                    {
                        from_exterior = !from_exterior;
                        ++isect_count;
                        if (bnd_extd.contains(point))
                        {
                            output.boundary_points.push_back(point);
                            output.directions.push_back(idir);
                            output.signs.push_back(sign_local);
                            sign_local = -sign_local;
                        }
                    };
                    
                    // Perform the ray trace along the [idir] axis in the [sign] direction
                    geom.trace_aligned_ray(idir, sign, x_comp, on_intersection);

                    
                    // if ((isect_count != 2) && (isect_count != 0) && (idir == 2))
                    // if ((isect_count != 2) && (isect_count != 0))
                    // {
                    //     print("DANG", isect_count);
                    //     print("g_t_i", global::tr_id);
                    //     auto& v = output.boundary_points.data(device::cpu);
                        
                    //     std::vector<pnt_t> pts;
                    //     pts.push_back(x_comp);
                    //     const auto on_intersection2 = [&](const auto& point)
                    //     {
                    //         pts.push_back(point);
                    //     };
                        
                    //     // Sad debugging
                    //     const auto& bvhh = geom.axis_bvh[idir];
                    //     geom.trace_aligned_ray(idir, sign, x_comp, on_intersection2);
                    //     print("npt incl. orig.:", pts.size());
                    //     io::output_vtk("bad_pts.vtk", pts, false);
                    //     geom::detail::debug_output_bvh("bvh_proj.vtk", bvhh);
                        
                    //     // output the triangles considered
                    //     using p2d_t = ctrs::array<real_t, 2>;
                    //     int t0 = (idir+1)%3;
                    //     int t1 = (idir+2)%3;
                    //     p2d_t x_prj = {x_comp[t0], x_comp[t1]};
                    //     std::vector<std::size_t> tris;
                    //     bvhh.check_elements([&](const auto& i) {tris.push_back(i);}, {x_comp[t0], x_comp[t1]});
                    //     geom::detail::output_subset("tris.vtk", geom, tris);
                    //     std::cin.get();
                    // }
                    
                    // global::tr_id++;
                };
                
                // Port to GPU later -_-
                dispatch::execute(ibound, loop_load, device::cpu);
            }
        }
        
        return output;
    }
}