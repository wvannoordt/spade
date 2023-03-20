#pragma once

#include <iostream>

#include "core/ctrs.h"

namespace spade::coords
{
    template <typename float_t> struct point_t : public ctrs::arithmetic_array_t<float_t, 3, point_t<float_t>>
    {
        float_t& x() { return (*this)[0]; }
        float_t& y() { return (*this)[1]; }
        float_t& z() { return (*this)[2]; }
        const float_t& x() const { return (*this)[0]; }
        const float_t& y() const { return (*this)[1]; }
        const float_t& z() const { return (*this)[2]; }
    };

    template <typename float_t> std::ostream& operator << (std::ostream& os, const point_t<float_t>& point)
    {
        os << "point: {x: " << point.x() << ", y: " << point.y() << ", z: " << point.z() << "}";
        return os;
    }
}