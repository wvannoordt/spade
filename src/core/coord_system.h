#pragma once

#include <type_traits>
#include <functional>
#include <concepts>
#include "core/ctrs.h"

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
        t.bad();
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
        dtype eta0, eta1, deta, f0, norm, alpha0, alpha1, beta0, beta1;
    };
}