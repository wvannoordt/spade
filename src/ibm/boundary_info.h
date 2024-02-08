#pragma once

#include "core/ctrs.h"
#include "grid/grid.h"

namespace spade::ibm
{
    template <const int num_layers, typename float_t, template <typename> typename container_t = device::basic_shared_vector>
    struct ghost_info_t
    {
        using idx_t            = grid::cell_idx_t;
        using pnt_t            = coords::point_t<float_t>;
        using vec_t            = ctrs::array<float_t, pnt_t::size()>;
        using image_type       = ghost_info_t<num_layers, float_t, utils::vec_image_t>;
        using const_image_type = ghost_info_t<num_layers, float_t, utils::const_vec_image_t>;
        
        using sbool = utils::sbool;
        
        template <typename data_t> using layer_arr_t    = container_t<ctrs::array<data_t, num_layers>>;
        template <typename data_t> using no_layer_arr_t = container_t<data_t>;
        
        // Note that irregular cells are implicitly calculable from this information.
        layer_arr_t   <idx_t> indices;          // i, j, k, lb of the ghost points
        no_layer_arr_t<pnt_t> boundary_points;  // Point at which the ray intersects the boundary
        no_layer_arr_t<vec_t> boundary_normals; // Normal vector at the boundary point
        layer_arr_t   <pnt_t> closest_points;   // Closest point on the boundary
        layer_arr_t   <vec_t> closest_normals;  // Normal vector at the closest point on the boundary
        no_layer_arr_t<int>   signs;            // sign of 1 means x[dir]_ghost < x[dir]_irreg
        layer_arr_t   <sbool> can_fill;         // indicates whether or not the ghost cell can be inserted into the domain.
        
        void transfer()
        {
            indices.transfer();
            boundary_points.transfer();
            boundary_normals.transfer();
            closest_points.transfer();
            closest_normals.transfer();
            signs.transfer();
            can_fill.transfer();
        }
    };
    
    template <typename float_t, template <typename> typename container_t = device::basic_shared_vector>
    struct diag_ghost_info_t
    {
        using idx_t            = grid::cell_idx_t;
        using pnt_t            = coords::point_t<float_t>;
        using vec_t            = ctrs::array<float_t, pnt_t::size()>;
        using image_type       = diag_ghost_info_t<float_t, utils::vec_image_t>;
        using const_image_type = diag_ghost_info_t<float_t, utils::const_vec_image_t>;
        
        using sbool = utils::sbool;
        
        // Note that irregular cells are implicitly calculable from this information.
        container_t<idx_t> indices;          // i, j, k, lb of the ghost points
        container_t<pnt_t> closest_points;   // Closest point on the boundary
        container_t<vec_t> closest_normals;  // Normal vector at the closest point on the boundary
        container_t<sbool> can_fill;         // indicates whether or not the ghost cell can be inserted into the domain.
        
        void transfer()
        {
            indices.transfer();
            closest_points.transfer();
            closest_normals.transfer();
            can_fill.transfer();
        }
    };
    
    
    template <
        const int grid_dim,
        const int nlayers,
        typename float_t,
        template <typename> typename container_t = device::basic_shared_vector
        >
    struct boundary_info_t
    {
        using pnt_t = coords::point_t<float_t>;
        using vec_t = ctrs::array<float_t, pnt_t::size()>;
        
        constexpr static int num_layers() { return nlayers; }
        constexpr static int dim()        { return grid_dim; }
        ctrs::array<ghost_info_t<nlayers, float_t, container_t>, grid_dim> aligned;
        diag_ghost_info_t<float_t, container_t> diags;
        
        void transfer()
        {
            for (auto& l: aligned) l.transfer();
            diags.transfer();
        }
    };
    
    template <typename array_t, typename geom_t, typename kernel_t>
    inline auto
    compute_boundary_info(const array_t& sampled_array, const geom_t& geom, const kernel_t&)
    {
        using omni_type = typename kernel_t::omni_type;
        static_assert(omni_type::template max_extent<0> == -omni_type::template min_extent<0>,
            "Stencil for boundary info must be symmetric");
        static_assert(omni_type::template min_extent<1> == 0,
            "Stencil for boundary info must be face-aligned (no tangential offsets)");
        static_assert(omni_type::template max_extent<1> == 0,
            "Stencil for boundary info must be face-aligned (no tangential offsets)");
        static_assert(omni_type::template min_extent<2> == 0,
            "Stencil for boundary info must be face-aligned (no tangential offsets)");
        static_assert(omni_type::template max_extent<2> == 0,
            "Stencil for boundary info must be face-aligned (no tangential offsets)");
        
        constexpr int num_layers = static_math::moddiv<omni_type::template max_extent<0>,2>::value + 1;
        
        using grid_t   = typename array_t::grid_type;
        using real_t   = typename grid_t::coord_type;
        using pnt_t    = typename grid_t::coord_point_type;
        using output_t = boundary_info_t<grid_t::dim(), num_layers, real_t>;
        
        using vec_t  = typename output_t::vec_t;
        
        static_assert(std::same_as<typename grid_t::coord_sys_type, coords::identity<typename grid_t::coord_type>>, "ghosts currently only working for identity coordinates");
        
        const auto& grid = sampled_array.get_grid();
        
        for (int d = 0; d < grid_t::dim(); ++d)
        {
            int required = utils::max(num_layers, 1);
            if (sampled_array.get_num_exchange(d) < required)
            {
                std::string msg = "Attempted to generate boundary info without enough exchange cells provided. (Required ";
                msg += std::to_string(num_layers);
                msg += ", but found ";
                msg += std::to_string(sampled_array.get_num_exchange(d));
                msg += " in direction ";
                msg += std::to_string(d);
                msg += ")";
                throw except::sp_exception(msg);
            }
        }
        
        // Select the blocks that intersect the boundary
        const auto is_intersect = [&](const auto& lb_loc)
        {
            auto bnd_with_exch = grid.get_bounding_box(lb_loc);
            const auto dx = grid.get_dx(lb_loc);
            for (int d = 0; d < grid.dim(); ++d)
            {
                // bounding box includes exchange cells
                bnd_with_exch.min(d) -= sampled_array.get_num_exchange(d)*dx[d];
                bnd_with_exch.max(d) += sampled_array.get_num_exchange(d)*dx[d];
            }
            return geom.box_contains_boundary(bnd_with_exch);
        };
        const auto lbs = grid.select_blocks(is_intersect, partition::local);
        
        output_t output;

        // We will only trace rays from the bounding box of the geometry
        const auto geom_bbx = geom.get_bounding_box().inflate(1.05);
        
        // We will exclude ghost cells outside the domain
        const auto domain_bbx = grid.get_bounds();
        
        // Looping over boxes with boundary intersections
        for (int idir = 0; idir < grid.dim(); ++idir)
        {
            int idir0 = (idir + 1) % grid.dim();
            int idir1 = (idir + 2) % grid.dim();
            // We will only modify the list in this direction
            auto& list = output.aligned[idir];
            
            const auto ng = sampled_array.get_num_exchange(idir);
            const auto nx = grid.get_num_cells(idir);
            
            for (const auto& lb: lbs)
            {
                // We will not consider xray intersection points outside of the domain
                const auto domain_boundary = grid.is_domain_boundary(lb);
                const auto dx = grid.get_dx(lb);
                
                // Note that there is some issue with intersection detection in exchange cells
                
                // Compute the indices that we will trace from
                bound_box_t<int, 3> ibound;
                for (int d = 0; d < grid.dim(); ++d)
                {
                    ibound.min(d) = -1;
                    ibound.max(d) = grid.get_num_cells(d) + 1;
                }
                ibound.min(idir) = 0;
                ibound.max(idir) = 1;

                const auto loop_load = [&](const auto& ii)
                {
                    grid::cell_idx_t icell_orig(ii[0], ii[1], ii[2], lb.value);
                    
                    // Note: this is performed in computational coordinates, so it is assumed that the geometry is in computational coordinates
                    auto x_comp  = grid.get_comp_coords(icell_orig);
                    x_comp[idir] = geom_bbx.max(idir);
                    
                    
                    // Now we trace the ray {x_comp, ray_vec}
                    
                    // We assume that we start out with a point outside
                    // the geometry (later modify this for internal flows)
                    
                    // We will do this at every intersection point
                    int isect_count = 0;
                    const int sign = -1;
                    const auto on_intersection = [&](const auto& point, const auto& normal)
                    {
                        ++isect_count;
                        
                        
                        const auto  gp_sign  = -utils::sign(normal[idir]*sign);
                        const auto  lb       = utils::tag[partition::local](icell_orig.lb());
                        const auto& b_point  = point;
                        const auto& bbox     = grid.get_bounding_box(lb);
                        const auto  dx       = grid.get_dx(lb);
                            
                        // Note that we include xray points inside the exchange cells in the tracing direction
                        auto bnd_extd = grid.get_bounding_box(lb);
                        
                        if (!(domain_boundary.min(idir))) bnd_extd.min(idir) -= sampled_array.get_num_exchange(idir)*dx[idir];
                        else                              bnd_extd.min(idir) -= gp_sign*0.49999*dx[idir]; 
                        if (!(domain_boundary.max(idir))) bnd_extd.max(idir) += sampled_array.get_num_exchange(idir)*dx[idir];
                        else                              bnd_extd.max(idir) -= gp_sign*0.49999*dx[idir];
                        
                        bnd_extd.min(idir0) -= dx[idir0];
                        bnd_extd.max(idir0) += dx[idir0];
                        
                        bnd_extd.min(idir1) -= dx[idir1];
                        bnd_extd.max(idir1) += dx[idir1];
                            
                        auto icell = icell_orig;
                        // Note that only the index value in the traced direction is wrong
                        icell[idir]  = (b_point[idir]-bbox.min(idir) + (ng+1)*dx[idir])/dx[idir];
                        icell[idir] -= (ng+1);
                        
                        const auto xc_comp = grid.get_comp_coords(icell);
                        const auto xdiff   = b_point[idir] - xc_comp[idir];
                        if (xdiff*gp_sign < 0.0) icell[idir] -= gp_sign;
                        
                        // At this point, icell contains the correct location of the ghost cell.
                        // We proceed to find the nearest point on the geometry
                        const auto x_ghost = grid.get_comp_coords(icell);
                        real_t search_radius = 1e-9;
                        for (const auto& dx_v: dx) search_radius = utils::max(search_radius, 2*ng*dx_v);
                        const auto nearest_boundary_point = geom.find_closest_boundary_point(x_ghost, search_radius);
                        
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
                            nv  = real_t(0.0)*nv;
                            nv += normal;
                        }
                        nv /= ctrs::array_norm(nv);
                        
                        if (bnd_extd.contains(point))
                        {
                            
                            list.signs.push_back(gp_sign);
                            list.indices.push_back(icell);
                            list.boundary_points.push_back(point);
                            list.boundary_normals.push_back(normal);
                            list.closest_normals.push_back(nv);
                            list.closest_points.push_back(nearest_boundary_point);
                            
                            auto& icells_recent   = list.indices.back();
                            auto& nvs_recent      = list.closest_normals.back();
                            auto& bndypts_recent  = list.closest_points.back();
                            
                            
                            for (int ilayer = 1; ilayer < num_layers; ++ilayer)
                            {
                                auto& icell_lyr  = icells_recent [ilayer];
                                auto& nv_lyr     = nvs_recent    [ilayer];
                                auto& bndypt_lyr = bndypts_recent[ilayer];
                                
                                const auto truncate = [](const int val, const int lo, const int hi)
                                {
                                    return utils::min(utils::max(lo, val), hi);
                                };
                                
                                icell_lyr.i(idir) -= gp_sign*ilayer;
                                icell_lyr.i(idir)  = truncate(icell_lyr.i(idir), -ng, nx+ng-1);
                                
                                const auto xc_comp_lyr = grid.get_comp_coords(icell_lyr);
                                real_t lyr_search_radius = 1e-9;
                                for (const auto& dx_v: dx) lyr_search_radius = utils::max(lyr_search_radius, (2+ilayer)*ng*dx_v);
                                
                                const auto nearest_boundary_point_lyr = geom.find_closest_boundary_point(xc_comp_lyr, lyr_search_radius);
                                
                                vec_t nv_lyr_bndy = 0.0;
                                nv_lyr_bndy += (nearest_boundary_point_lyr - xc_comp_lyr);
                                
                                auto dist_lyr = ctrs::array_norm(nv_lyr_bndy);
                                if (dist_lyr < tol*diag)
                                {
                                    nv_lyr_bndy = real_t(0.0)*nv_lyr_bndy;
                                    nv_lyr_bndy += nv_lyr[0];
                                }
                                
                                nv_lyr_bndy /= ctrs::array_norm(nv_lyr_bndy);
                                bndypt_lyr = nearest_boundary_point_lyr;
                            }
                            
                            list.can_fill.push_back(utils::sbool{true});
                            auto& can_fill_recent = list.can_fill.back();
                        }
                    };
                    
                    // Perform the ray trace along the [idir] axis in the [sign] direction
                    geom.trace_aligned_ray(idir, sign, x_comp, on_intersection);
                };
                // Port to GPU later -_-
                dispatch::execute(ibound, loop_load, device::cpu);
            }
            
            std::size_t num_in_this_dir = list.indices.size();
            for (std::size_t id = 0; id < num_in_this_dir; ++id)
            {
                for (int ilayer = 0; ilayer < num_layers; ++ilayer)
                {
                    const auto icell_loc = list.indices[id][ilayer];
                    const auto xg        = grid.get_comp_coords(icell_loc);
                    int i_idir = icell_loc.i(idir);
                    const auto domain_boundary = grid.is_domain_boundary(utils::tag[partition::local](icell_loc.lb()));
                    bool is_domain_bndy_cell = (domain_boundary.min(idir) && (i_idir < 0)) || (domain_boundary.max(idir) && (i_idir >= nx));
                    list.can_fill[id][ilayer]  = geom.is_interior(xg);// || is_domain_bndy_cell;
                }
            }
        }
        output.transfer();
        return output;
    }
}