#pragma once

#include <concepts>
#include "ctrs.h"

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
    
    template <typename dtype> struct identity
    {
        typedef dtype coord_type;
        template <ctrs::vec_nd<3,coord_type> coords_t> auto map(const coords_t& coords) const {return coords;}
    };
}