#pragma once

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "core/except.h"
#include "core/ctrs.h"
#include "core/point.h"
#include "core/ascii.h"
#include "core/except.h"

#include "io/io_vtk.h"

#include "geom/primitives.h"
#include "geom/bvh.h"

namespace spade::geom
{
    template <const int dim, const int n_edge = 3, typename float_t = double>
    struct vtk_geom_t
    {
        using value_type = float_t;
        using uint_t     = std::size_t;
        template <typename T> using T_pnt_t = coords::point_t<T>;
        using pnt_t      = T_pnt_t<float_t>;
        using vec_t      = ctrs::array<float_t,3>;
        using face_t     = ctrs::array<uint_t, n_edge>;
        
        static constexpr int num_edges() { return n_edge; }
        
        bool                       is_external; //true -> points at infinity are on the exterior
        std::vector<pnt_t>         points;
        std::vector<vec_t>         normals; // These will ALWAYS point from the boundary towards the computational/interior domain
        std::vector<face_t>        faces;
        bound_box_t<float_t,dim>   bbox, bbox_inflated;
        bvh_impl_t<3, float_t, std::vector> vol_bvh;
        bvh_impl_t<2, float_t, std::vector> axis_bvh[3];
        std::size_t max_2d_bvh_elems;
        
        pnt_t centroid(const uint_t& idx) const
        {
            pnt_t output;
            const auto& face = faces[idx];
            const auto& p0   = points[face[0]];
            const auto& p1   = points[face[1]];
            const auto& p2   = points[face[2]];
            
            #pragma unroll
            for (int d = 0; d < output.size(); ++d)
            {
                output[d] = (p0[d] + p1[d] + p2[d])/3.0;
            }
            return output;
        }
        
        auto get_bounding_box() const
        {
            return bbox;
        }
        
        void calc_bvhs()
        {
            const auto box_check = [&](const uint_t& i, const auto& bnd)
            {
                const auto& face = faces[i];
                return primitives::box_tri_intersect(points[face[0]], points[face[1]], points[face[2]], bnd);
            };
            
            vol_bvh.calculate(bbox_inflated, faces.size(), box_check);
            
            max_2d_bvh_elems = 0;
            for (int d = 0; d < 3; ++d)
            {
                auto& bvh    = axis_bvh[d];
                const int d0 = (d+1)%3;
                const int d1 = (d+2)%3;
                
                const auto tol = 1e-5;
                const auto box_check = [&](const uint_t& i, const auto& bnd)
                {
                    bound_box_t<float_t, 3> bbox_3d;
                    bbox_3d.min(d0)    = bnd.min(0);
                    bbox_3d.max(d0)    = bnd.max(0);
                    bbox_3d.min(d1)    = bnd.min(1);
                    bbox_3d.max(d1)    = bnd.max(1);
                    bbox_3d.min(d)     = bbox_inflated.min(d);
                    bbox_3d.max(d)     = bbox_inflated.max(d);
                    const auto& face   = faces[i];
                    const auto& normal = normals[i];
                    auto plane_normal  = normals[i]; //avoid decltype shenanigans
                    plane_normal = 0.0;
                    plane_normal[d] = 1.0;
                    if (utils::abs(ctrs::dot_prod(normal, plane_normal))<tol) return false;
                    return primitives::box_tri_intersect(points[face[0]], points[face[1]], points[face[2]], bbox_3d);
                };
                
                bound_box_t<float_t, 2> bbox_2d;
                bbox_2d.min(0) = bbox_inflated.min(d0);
                bbox_2d.max(0) = bbox_inflated.max(d0);
                bbox_2d.min(1) = bbox_inflated.min(d1);
                bbox_2d.max(1) = bbox_inflated.max(d1);
                
                bvh.calculate(bbox_2d, faces.size(), box_check);
                max_2d_bvh_elems = utils::max(max_2d_bvh_elems, bvh.get_max_elems());
            }
            
        }
        
        template <typename pt_flt_t>
        T_pnt_t<pt_flt_t> find_closest_boundary_point(const T_pnt_t<pt_flt_t>& x_in, const pt_flt_t& search_radius) const
        {
            // if (!(bbox_inflated.contains(x_in)))
            // {
            //     throw except::sp_exception("can't find closest point if outside bvh!");
            // }
            
            const int bvhdim  = decltype(vol_bvh)::bvh_dim();
            bound_box_t<float_t, bvhdim> search_box;
            for (int i = 0; i < vol_bvh.bvh_dim(); ++i)
            {
                search_box.min(i) = x_in[i] - search_radius;
                search_box.max(i) = x_in[i] + search_radius;
            }
            
            
            pt_flt_t min_dist = 1e50;
            T_pnt_t<pt_flt_t> output(-1e9,-1e9,-1e9);
            const auto point_check = [&](const auto& iface)
            {
                const auto& face  = faces[iface];
                const auto& p0    = points[face[0]];
                const auto& p1    = points[face[1]];
                const auto& p2    = points[face[2]];
                const auto  x_tri = primitives::closest_point_on_tri(x_in, p0, p1, p2);
                const auto  dist  = ctrs::array_norm(x_tri - x_in);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    output = x_tri;
                }
            };
            
            vol_bvh.check_elements(point_check, search_box);            
            return output;
        }
        
        template <typename rhs_float_t>
        bool box_contains_boundary(const bound_box_t<rhs_float_t, 3>& bnd) const
        {
            bool output = false;
            const auto check_tri = [&](const auto& i)
            {
                const auto& face   = faces[i];
                const bool success = primitives::box_tri_intersect(points[face[0]], points[face[1]], points[face[2]], bnd);
                return success;
            };

            vol_bvh.check_elements([&](const auto& i) { output = output || check_tri(i); }, bnd);

            return output;
        }
        
        template <const int ar_size, const int rec_size, typename pfloat_t, typename on_isect_t>
        void trace_aligned_ray_rec_impl(const int dir, const int sign, const coords::point_t<pfloat_t>& x, const on_isect_t& on_isect) const
        {
            using vec_t = ctrs::array<pfloat_t, 3>;
            using p2d_t = ctrs::array<pfloat_t, 2>;
            vec_t nvec = 0.0;
            nvec[dir] = sign;
            
            int t0 = (dir+1)%3;
            int t1 = (dir+2)%3;
            
            p2d_t x_prj = {x[t0], x[t1]};
            
            const auto& table = axis_bvh[dir];
            using bvh_type = utils::remove_all<decltype(table)>::type;
            using bvh_pt_t = typename bvh_type::pnt_t;
            bvh_pt_t point_2d;
            point_2d[0] = x[t0];
            point_2d[1] = x[t1];
            
            ctrs::array<pnt_t,  ar_size> xray_buf;
            ctrs::array<uint_t, ar_size> tri_buf;
            ctrs::array<uint_t, ar_size> permute_buf;
            
            uint_t buf_idx = 0;
            bool overflow = false;
            const auto push_back = [&](const pnt_t& xx, const uint_t trii, const uint_t perm)
            {
                xray_buf[buf_idx]    = xx;
                tri_buf[buf_idx]     = trii;
                permute_buf[buf_idx] = perm;
                ++buf_idx;
                overflow = (buf_idx == ar_size);
                buf_idx = utils::min(buf_idx, ar_size - 1);
            };
            
            int isct = 0;
            const auto isect_check = [&](const auto& i)
            {
                ++isct;
                const auto& face = faces  [i];
                const auto& norm = normals[i];
                
                const auto& p0 = points[face[0]];
                const auto& p1 = points[face[1]];
                const auto& p2 = points[face[2]];
                
                p2d_t p0_prj = {p0[t0], p0[t1]};
                p2d_t p1_prj = {p1[t0], p1[t1]};
                p2d_t p2_prj = {p2[t0], p2[t1]};
                
                if (primitives::point_in_tri(x_prj, p0_prj, p1_prj, p2_prj))
                {
                    const auto t = sign*ctrs::dot_prod(x-p0, norm)/norm[dir];
                    pnt_t xb = x;
                    xb -= t*nvec;
                    push_back(xb, i, buf_idx);
                }
            };
            table.check_elements(isect_check, point_2d);
            if (overflow)
            {
                if constexpr (rec_size > 0) trace_aligned_ray_rec_impl<2*ar_size, rec_size-1>(dir, sign, x, on_isect);
                else
                {
                    throw except::sp_exception("Too many xray points traced! Found more than " + std::to_string(ar_size) + " intersections");
                }
            }
            else
            {
                if (buf_idx == 0) return;
                std::sort(permute_buf.begin(), permute_buf.begin() + buf_idx, [&](const auto& a, const auto& b){ return sign*xray_buf[a][dir] < sign*xray_buf[b][dir]; });
                int tsgn = utils::sign(normals[tri_buf[permute_buf[0]]][dir]);
                on_isect(xray_buf[permute_buf[0]], normals[tri_buf[permute_buf[0]]]);
                for (int j = 1; j < buf_idx; ++j)
                {
                    if (tsgn*utils::sign(normals[tri_buf[permute_buf[j]]][dir]) < 0.0)
                    {
                        tsgn = -tsgn;
                        on_isect(xray_buf[permute_buf[j]], normals[tri_buf[permute_buf[j]]]);
                    }
                }
            }
        }
        
        template <typename pfloat_t, typename on_isect_t>
        void trace_aligned_ray(const int dir, const int sign, const coords::point_t<pfloat_t>& x, const on_isect_t& on_isect) const
        {
            trace_aligned_ray_rec_impl<128, 6>(dir, sign, x, on_isect);
        }
        
        bool is_interior(const pnt_t& x) const
        {
            int rcount = 0;
            constexpr int dir = 0;
            this->trace_aligned_ray(dir, 1, x, [&](const auto& xxx, const auto&)
            {
                if (xxx[dir] > x[dir]) ++rcount;
            });
            return (rcount % 2) != 0;
        }
    };
    
    template <const int dim, const int n_edge, typename float_t>
    static std::string read_vtk_geom(const std::string& fname, vtk_geom_t<dim, n_edge, float_t>& surf, const bool is_external)
    {
        using uint_t = typename vtk_geom_t<dim, n_edge, float_t>::uint_t;
        using face_t = typename vtk_geom_t<dim, n_edge, float_t>::face_t;
        using pnt_t  = typename vtk_geom_t<dim, n_edge, float_t>::pnt_t;
        using vec_t  = typename vtk_geom_t<dim, n_edge, float_t>::vec_t;
        
        surf.is_external = is_external;
        
        utils::ascii_file_t fh(fname);
        fh.next_line();
        std::string header = fh.next_line();
        fh.next_line();
        fh.expect("ASCII");
        
        fh.next_line();
        fh.expect("DATASET POLYDATA");
        
        fh.next_line();
        std::string j0, j1;
        std::size_t num_points;
        fh.parse(j0, num_points, j1);
        
        surf.bbox.min(0) =  1e50;
        surf.bbox.min(1) =  1e50;
        surf.bbox.min(2) =  1e50;
        surf.bbox.max(0) = -1e50;
        surf.bbox.max(1) = -1e50;
        surf.bbox.max(2) = -1e50;
        surf.points.resize(num_points);
        std::size_t i_pt = 0;
        int i_xyz = 0;
        float_t x;
        while (i_pt < num_points)
        {
            while (!fh.try_parse(x)) fh.next_line();
            surf.points[i_pt][i_xyz] = x;
            i_xyz++;
            if (i_xyz >= 3)
            {
                i_xyz = 0;
                i_pt++;
            }
        }

        for (std::size_t i = 0; i < num_points; ++i)
        {
            const auto& pt = surf.points[i];

            surf.bbox.min(0) = utils::min(surf.bbox.min(0), pt[0]);
            surf.bbox.min(1) = utils::min(surf.bbox.min(1), pt[1]);
            surf.bbox.min(2) = utils::min(surf.bbox.min(2), pt[2]);
            surf.bbox.max(0) = utils::max(surf.bbox.max(0), pt[0]);
            surf.bbox.max(1) = utils::max(surf.bbox.max(1), pt[1]);
            surf.bbox.max(2) = utils::max(surf.bbox.max(2), pt[2]);
        }
        
        surf.bbox_inflated = surf.bbox.inflate(1.0001);
        
        std::size_t num_faces, junk;
        fh.next_line();
        fh.parse(j0, num_faces, junk);

        uint_t fsiz;
        face_t floc;
        surf.faces.reserve(num_faces);
        surf.normals.reserve(num_faces);
        for (std::size_t i = 0; i < num_faces; ++i)
        {
            fh.next_line();
            fh.parse(fsiz, floc);
            surf.faces.push_back(floc);
            pnt_t pp0 = surf.points[floc[0]];
            pnt_t pp1 = surf.points[floc[1]];
            pnt_t pp2 = surf.points[floc[2]];
            vec_t p0;
            vec_t p1;
            vec_t p2;
            ctrs::copy_array(pp0, p0);
            ctrs::copy_array(pp1, p1);
            ctrs::copy_array(pp2, p2);
            vec_t n = ctrs::cross_prod(p1-p0, p2-p1);
            n /= ctrs::array_norm(n);
            surf.normals.push_back(n);
        }
        
        surf.calc_bvhs();
        
        return header;
    }
    
    namespace detail
    {
        template <typename vtk_t, typename subset_t>
        static void output_subset(const std::string& fname, const vtk_t& geom, const subset_t& subset)
        {
            std::ofstream mf(fname);
            std::size_t count = 0;
            for (const auto& ifc: subset) count++;
            
            mf << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << 3*count << " double\n";
            for (const auto& ifc: subset)
            {
                const auto& fc = geom.faces[ifc];
                const auto& p0 = geom.points[fc[0]];
                const auto& p1 = geom.points[fc[1]];
                const auto& p2 = geom.points[fc[2]];
                mf << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
                mf << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
                mf << p2[0] << " " << p2[1] << " " << p2[2] << "\n";
            }
            mf << "POLYGONS " << count << " " << 4*count << "\n";
            std::size_t idx = 0;
            for (std::size_t i = 0; i < count; ++i)
            {
                mf << "3 " << idx++ << " " << idx++ << " " << idx++ << "\n";
            }
            
            mf << "CELL_DATA " << count << "\n";
            mf << "SCALARS Data double\nLOOKUP_TABLE default\n";
            for (std::size_t i = 0; i < count; ++i)
            {
                mf << "1" << "\n";
            }
            
        }
    }
}