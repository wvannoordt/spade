#pragma once

#include "core/ctrs.h"
#include "core/bounding_box.h"
#include "core/utils.h"

namespace spade::geom::primitives
{
    template <typename T> concept point_3d = (T::size() == 3 && std::floating_point<typename T::value_type>);
    template <typename T> concept point_2d = (T::size() == 2 && std::floating_point<typename T::value_type>);
    
    template <point_3d arr_t, point_3d nv_t>
    _sp_hybrid static bool box_plane_intersect(const arr_t& point, const nv_t& normalVec, const bound_box_t<typename arr_t::value_type, 3>& box)
    {
        //strictly speaking, only need to test the "most aligned" vertices
        using float_t = typename arr_t::value_type;
        bool hasPositive = false;
	    bool hasNegative = false;
        float_t corner[3];
        float_t dotProd;
        for (char b = 0; b < 8; b++)
        {
            corner[0] = box(0,((b>>0)&1));
            corner[1] = box(1,((b>>1)&1));
            corner[2] = box(2,((b>>2)&1));
            dotProd  = normalVec[0]*(corner[0] - point[0]);
            dotProd += normalVec[1]*(corner[1] - point[1]);
            dotProd += normalVec[2]*(corner[2] - point[2]);
            hasPositive = hasPositive||(dotProd>0);
            hasNegative = hasNegative||(dotProd<0);
            if ((float_t(-1e-8) < dotProd) && (dotProd < float_t(1e-8))) return true;
            if (hasNegative&&hasPositive) return true;
        }
        return false;
    }
        
    
    template <point_3d arr_t>
    _sp_hybrid static bool box_tri_intersect(
        const arr_t& p1,
        const arr_t& p2,
        const arr_t& p3,
        const bound_box_t<typename arr_t::value_type, 3>& box)
    {
        using v_t     = arr_t;
        using float_t = typename arr_t::value_type;
        using bnd_t   = bound_box_t<float_t, 3>;
        float_t half{0.5};
        v_t c{half*(box.min(0) + box.max(0)), half*(box.min(1) + box.max(1)), half*(box.min(2) + box.max(2))};
        v_t h{half*box.size(0), half*box.size(1), half*box.size(2)};
        bnd_t tbox;
        for (int i = 0; i < v_t::size(); ++i)
        {
            tbox.min(i) = utils::min(p1[i], p2[i], p3[i]);
            tbox.max(i) = utils::max(p1[i], p2[i], p3[i]);
        }
        
        //TODO v  v  v  v
        auto nvec = ctrs::cross_prod(p2-p1, p3-p2);
        if (!tbox.intersects(box)) return false;
        if (!box_plane_intersect(p1, nvec, box)) return false;
        float_t ei[3];
        float_t fj[3][3];
        fj[0][0] = p2[0] - p1[0];
        fj[1][0] = p2[1] - p1[1];
        fj[2][0] = p2[2] - p1[2];
        fj[0][1] = p3[0] - p2[0];
        fj[1][1] = p3[1] - p2[1];
        fj[2][1] = p3[2] - p2[2];
        fj[0][2] = p1[0] - p3[0];
        fj[1][2] = p1[1] - p3[1];
        fj[2][2] = p1[2] - p3[2];
        float_t aij[3];
        float_t p[3];
        float_t r;
        int i, j;
        float_t pmin, pmax;
        for (int z = 0; z < 9; z++)
        {
            i = z/3;
            j = z%3;
            ei[0] = 0.0; ei[1] = 0.0; ei[2] = 0.0;
            ei[i] = 1.0;
            aij[0] = ei[1]*fj[2][j] - ei[2]*fj[1][j];
            aij[1] = ei[2]*fj[0][j] - ei[0]*fj[2][j];
            aij[2] = ei[0]*fj[1][j] - ei[1]*fj[0][j];
            p[0] = aij[0]*(p1[0]-c[0]) + aij[1]*(p1[1]-c[1]) + aij[2]*(p1[2]-c[2]);
            p[1] = aij[0]*(p2[0]-c[0]) + aij[1]*(p2[1]-c[1]) + aij[2]*(p2[2]-c[2]);
            p[2] = aij[0]*(p3[0]-c[0]) + aij[1]*(p3[1]-c[1]) + aij[2]*(p3[2]-c[2]);
            pmin = utils::min(p[0], p[1], p[2]);
            pmax = utils::max(p[0], p[1], p[2]);
            r = h[0]*utils::abs(aij[0]) + h[1]*utils::abs(aij[1]) + h[2]*utils::abs(aij[2]);
            if ((r < pmin) || (pmax < -r)) return false;
        }
        return true;
    }
    
    template <point_2d arr_t>
    _sp_hybrid static bool point_in_tri(const arr_t& x, const arr_t& p0, const arr_t& p1, const arr_t& p2, const typename arr_t::value_type tol = 1e-7)
    {
        using real_t = typename arr_t::value_type;
        const auto tri_area = [](const arr_t& a, const arr_t& b, const arr_t& c)
        {
            auto dx0 = b[0] - a[0];
            auto dy0 = b[1] - a[1];
            auto dx1 = c[0] - b[0];
            auto dy1 = c[1] - b[1];
            
            // i    j    k
            // dx0  dy0  0
            // dx1  dy1  0
            
            return dx0*dy1 - dy0*dx1;
        };
        
        const auto a_all = utils::abs(tri_area(p0, p1, p2));
        const auto a_0   = utils::abs(tri_area(x,  p1, p2));
        const auto a_1   = utils::abs(tri_area(p0, x,  p2));
        const auto a_2   = utils::abs(tri_area(p0, p1, x ));
        
        return utils::abs(a_all - a_0 - a_1 - a_2) < tol;
    }
}