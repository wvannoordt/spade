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

#include "geom/primitives.h"
#include "geom/bvh.h"

namespace spade::geom
{
    template <const int dim, const int n_edge = 3, typename float_t = double>
    struct vtk_geom_t
    {
        using value_type = float_t;
        using uint_t     = std::size_t;
        using pnt_t      = coords::point_t<float_t>;
        using vec_t      = ctrs::array<float_t,3>;
        using face_t     = ctrs::array<uint_t, n_edge>;
        
        static constexpr int num_edges() { return n_edge; }
        
        std::vector<pnt_t>         points;
        std::vector<vec_t>         normals;
        std::vector<face_t>        faces;
        bound_box_t<float_t,dim>   bbox, bbox_inflated;
        
        
        bvh_t<3, float_t, std::vector> vol_bvh;
        
        void calc_bvh()
        {
            const auto box_check = [&](const uint_t& i, const auto& bnd)
            {
                const auto& face = faces[i];
                return primitives::box_tri_intersect(points[face[0]], points[face[1]], points[face[2]], bnd);
            };
            
            const auto tri_bbx = [&](const uint_t& i)
            {
                const auto& face = faces[i];
                bound_box_t<float_t,dim> output;
                
                for (int d = 0; d < dim; ++d)
                {
                    output.min(d) =  1e50;
                    output.max(d) = -1e50;
                }
                
                for (int k = 0; k < face_t::size(); ++k)
                {
                    const auto& pt = points[face[k]];
                    for (int d = 0; d < dim; ++d)
                    {                    
                        output.min(d) = utils::min(output.min(d), pt[d]);
                        output.max(d) = utils::max(output.max(d), pt[d]);
                    }
                }
            };
            vol_bvh.calculate(bbox_inflated, faces.size(), box_check, tri_bbx);
            detail::debug_output_bvh("bvh.vtk", vol_bvh);
        }
        
        template <typename rhs_float_t>
        bool partially_contained_by(const bound_box_t<rhs_float_t, 3>& bnd) const
        {
            //This is a total garbage implementation but I just can't be bothered to improve it at the moment
            for (std::size_t i = 0; i < faces.size(); ++i)
            {
                const auto& face   = faces[i];
                const bool success = primitives::box_tri_intersect(points[face[0]], points[face[1]], points[face[2]], bnd);
                if (success) return true;
            }
            return false;
        }
    };
    
    template <const int dim, const int n_edge, typename float_t>
    static std::string read_vtk_geom(const std::string& fname, vtk_geom_t<dim, n_edge, float_t>& surf)
    {
        using uint_t = typename vtk_geom_t<dim, n_edge, float_t>::uint_t;
        using face_t = typename vtk_geom_t<dim, n_edge, float_t>::face_t;
        using pnt_t  = typename vtk_geom_t<dim, n_edge, float_t>::pnt_t;
        using vec_t  = typename vtk_geom_t<dim, n_edge, float_t>::vec_t;
        
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
        
        float_t x, y, z;
        surf.bbox.min(0) =  1e50;
        surf.bbox.min(1) =  1e50;
        surf.bbox.min(2) =  1e50;
        surf.bbox.max(0) = -1e50;
        surf.bbox.max(1) = -1e50;
        surf.bbox.max(2) = -1e50;
        surf.points.reserve(num_points);
        for (std::size_t i = 0; i < num_points; ++i)
        {
            fh.next_line();
            fh.parse(x, y, z);
            surf.points.push_back(pnt_t(x, y, z));
            surf.bbox.min(0) = utils::min(surf.bbox.min(0), x);
            surf.bbox.min(1) = utils::min(surf.bbox.min(1), y);
            surf.bbox.min(2) = utils::min(surf.bbox.min(2), z);
            surf.bbox.max(0) = utils::max(surf.bbox.max(0), x);
            surf.bbox.max(1) = utils::max(surf.bbox.max(1), y);
            surf.bbox.max(2) = utils::max(surf.bbox.max(2), z);
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
        
        surf.calc_bvh();
        
        return header;
    }
}