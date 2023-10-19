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
        using vec_t  = typename ghost_list_t<typename grid_t::coord_type>::vec_t;
        
        
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
                
                // Note that there is some issue with intersection detection in exchange cells
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
                    const auto on_intersection = [&](const auto& point, const auto& normal)
                    {
                        from_exterior = !from_exterior;
                        sign_local    = -sign_local;
                        ++isect_count;
                        
                        if (bnd_extd.contains(point))
                        {
                            output.boundary_points.push_back(point);
                            output.directions.push_back(idir);
                            output.signs.push_back(-utils::sign(normal[idir]*sign));
                            output.indices.push_back(icell);
                            output.boundary_normals.push_back(normal);
                        }
                    };
                    
                    // Perform the ray trace along the [idir] axis in the [sign] direction
                    geom.trace_aligned_ray(idir, sign, x_comp, on_intersection);
                };
                
                
                // Port to GPU later -_-
                dispatch::execute(ibound, loop_load, device::cpu);
            }
        }
        
        // By this point, we have performed all of the x-ray tracing to get the surface points.
        // now, we will go and figure out the ghost cells from each boundary point.
        
        // Note that the indices in output.indices currently stores the cell that 'shot' the ray
        // so only the block is correct. We will use this to compute the correct index
        
        const auto& v = output.boundary_points.data(device::cpu);
        dispatch::ranges::linear_range_t rg(0UL, v.size(), device::cpu);
        const auto comp_ghost_info = [&](const auto& i)
        {
            auto& icell         = output.indices[i];
            const auto  lb      = utils::tag[partition::local](icell.lb());
            const auto& b_point = output.boundary_points[i];
            const auto& bbox    = grid.get_bounding_box(lb);
            const auto& idir    = output.directions[i];
            const auto& sign    = output.signs[i];
            const auto  dx      = grid.get_dx(lb);
            const auto  ng      = grid.get_num_exchange(idir);
            
            // Note that only the index value in the traced direction is wrong
            icell[idir]  = (b_point[idir]-bbox.min(idir) + (ng+1)*dx[idir])/dx[idir];
            icell[idir] -= (ng+1);
                        
            const auto xc_comp = grid.get_comp_coords(icell);
            const auto xdiff   = b_point[idir] - xc_comp[idir];
            if (xdiff*sign < 0.0) icell[idir] -= sign;
            
            // At this point, icell contains the correct location of the ghost cell.
            // We procede to find the nearest point on the geometry
            const auto x_ghost = grid.get_comp_coords(icell);
            real_t search_radius = 1e-9;
            for (const auto& dx_v: dx) search_radius = utils::max(search_radius, dx_v);
            const auto nearest_boundary_point = geom.find_closest_boundary_point(x_ghost, search_radius);
            output.closest_points.push_back(nearest_boundary_point);
            
            // Compute direction to the closest boundary point
            vec_t nv = 0.0;
            nv += (nearest_boundary_point - x_ghost);
            auto dist       = ctrs::array_norm(nv);
            const auto diag = ctrs::array_norm(dx);
            
            // Figure out this magic number
            const auto tol = 5e-3;
            
            // If the ghost point is too close to the boundary, then
            // the normal vector is just chosen as the surface normal vector
            if (dist < tol*diag)
            {
                nv  = 0.0*nv;
                nv += output.boundary_normals[i];
            }
            nv /= ctrs::array_norm(nv);
            output.closest_normals.push_back(nv);
        };
        
        dispatch::execute(rg, comp_ghost_info, device::cpu);
        
        output.transfer();
        
        return output;
    }
}