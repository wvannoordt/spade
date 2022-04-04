#pragma once

#include <type_traits>
#include <functional>
#include <concepts>
#include "core/ctrs.h"
#include "core/grid_index_types.h"

namespace cvdf::coords
{
    template <class T> concept coordinate_system = std::floating_point<typename T::coord_type>
    && requires(T t, ctrs::array<typename T::coord_type,3> x)
    {
        t;
        { t.map(x) } -> ctrs::vec_nd<3,typename T::coord_type>;
    };
    
    template <class T> concept diagonal_coordinate_system = coordinate_system<T> && requires(T t)
    {
        t;
        t.xcoord;
        t.ycoord;
        t.zcoord;
    };
    
    template <class T> concept coord_mapping_1D = requires(T t, typename T::coord_type dat)
    {
        dat = t.map(dat);
    };
    
    template <class T> concept coord_mapping_callable_1D = requires(T t)
    {
        t(0.0);
    };
    
    template <typename dtype> struct identity
    {
        typedef dtype coord_type;
        template <ctrs::vec_nd<3,coord_type> coords_t> auto map(const coords_t& coords) const {return coords;}
    };
    
    template <typename dtype> struct identity_1D
    {
        typedef dtype coord_type;
        dtype map(const dtype& coord) const {return coord;}
    };
    
    template <
        coord_mapping_1D xcoord_t,
        coord_mapping_1D ycoord_t,
        coord_mapping_1D zcoord_t>
    struct diagonal_coords
    {
        typedef decltype(typename xcoord_t::coord_type() + typename ycoord_t::coord_type() + typename zcoord_t::coord_type()) coord_type;
        diagonal_coords(void){}
        diagonal_coords(const xcoord_t& xcoord_in, const ycoord_t& ycoord_in, const zcoord_t& zcoord_in)
        {
            xcoord = xcoord_in;
            ycoord = ycoord_in;
            zcoord = zcoord_in;
        }
        template <ctrs::vec_nd<3,coord_type> coords_t> auto map(const coords_t& coords) const
        {
            return coords_t(xcoord.map(coords[0]), ycoord.map(coords[1]), zcoord.map(coords[2]));
        }
        xcoord_t xcoord;
        ycoord_t ycoord;
        zcoord_t zcoord;
    };
    
    template <coord_mapping_callable_1D callable_t>
    struct analytical_1D
    {
        typedef decltype(((callable_t*)(NULL))->operator()(0.0)) coord_type;
        typedef decltype(((callable_t*)(NULL))->operator()(0.0)) dtype;
        analytical_1D(){}
        analytical_1D(const callable_t& func_in) {func = &func_in;}
        dtype map(const dtype& coord) const {return (*func)(coord);}
        const callable_t* func;
    };
    
    template <typename dtype> struct integrated_tanh_1D
    {
        typedef dtype coord_type;
        integrated_tanh_1D(void){}
        integrated_tanh_1D(const dtype& y0, const dtype& y1, const dtype& inflation, const dtype& rate)
        {
            eta1 = y0;
            eta0 = y1;
            deta = eta1 - eta0;
            const dtype k = rate;
            const dtype strch = inflation;
            alpha0 = k/deta;
            alpha1 = -k/deta;
            beta0  =  alpha0*(strch*deta-eta0);
            beta1  =  alpha0*(eta1+strch*deta);
            auto abs = [](const dtype& d) -> dtype {return d<0?-d:d;};
            auto func = [&](const dtype& eta) -> dtype
            {
                return  log(abs(cosh(alpha0*eta+beta0)))/alpha0  +  log(abs(cosh(alpha1*eta+beta1)))/alpha1 - eta;
            };
            f0 = func(eta0);
            norm = func(eta1)-f0;
        }
        dtype map(const dtype& coord) const
        {
            auto abs = [](const dtype& d) -> dtype {return d<0?-d:d;};
            auto func = [&](const dtype& eta) -> dtype
            {
                return  log(abs(cosh(alpha0*eta+beta0)))/alpha0  +  log(abs(cosh(alpha1*eta+beta1)))/alpha1 - eta;
            };
            return eta0 + deta*(func(coord) - f0)/norm;
        }
        dtype metric(const dtype& coord) const
        {
            return norm/(deta*(tanh(alpha0*coord+beta0) + tanh(alpha1*coord+beta1) - 1.0));
        }
        dtype eta0, eta1, deta, f0, norm, alpha0, alpha1, beta0, beta1;
    };
    
    template <typename dtype, typename integral_t>
    ctrs::array<dtype, 3> calc_normal_vector(
        const identity<dtype>& coord,
        const ctrs::array<dtype, 3>& coords,
        const ctrs::array<grid::cell_t<integral_t>, 4>& i,
        const integral_t& idir)
    {
        ctrs::array<dtype, 3> output(0.0, 0.0, 0.0);
        output[idir]=1.0;
        return output;
    }
    
    template <typename dtype, typename integral_t>
    dtype calc_jacobian(
        const identity<dtype>& coord,
        const ctrs::array<dtype, 3>& coords,
        const ctrs::array<grid::cell_t<integral_t>, 4>& i)
    {
        return 1.0;
    }
    
    template <diagonal_coordinate_system coord_t, typename integral_t>
    ctrs::array<typename coord_t::coord_type, 3> calc_normal_vector(
        const coord_t& coord,
        const ctrs::array<typename coord_t::coord_type, 3>& coords,
        const ctrs::array<grid::cell_t<integral_t>, 4>& i,
        const integral_t& idir)
    {
        ctrs::array<typename coord_t::coord_type, 3> output(0.0, 0.0, 0.0);
        output[idir] = 1.0;
        
        output[0]*=coord.xcoord.metric(coords[0]);
        output[1]*=coord.ycoord.metric(coords[1]);
        output[2]*=coord.zcoord.metric(coords[2]);
        return output;
    }
    
    template <diagonal_coordinate_system coord_t, typename integral_t>
    typename coord_t::coord_type calc_jacobian(
        const coord_t& coord,
        const ctrs::array<typename coord_t::coord_type, 3>& coords,
        const ctrs::array<grid::cell_t<integral_t>, 4>& i)
    {
        return coord.xcoord.metric(coords[0])*coord.ycoord.metric(coords[1])*coord.zcoord.metric(coords[2]);
    }
}