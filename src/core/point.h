#pragma once

#include <iostream>

#include "core/ctrs.h"

namespace spade::coords
{
    template <typename float_t>
    struct point_t : public ctrs::arithmetic_array_t<float_t, 3, point_t<float_t>>
    {
        using base_t = ctrs::arithmetic_array_t<float_t, 3, point_t<float_t>>;
        using base_t::base_t;
        
        _sp_hybrid float_t& x() { return (*this)[0]; }
        _sp_hybrid float_t& y() { return (*this)[1]; }
        _sp_hybrid float_t& z() { return (*this)[2]; }
        _sp_hybrid const float_t& x() const { return (*this)[0]; }
        _sp_hybrid const float_t& y() const { return (*this)[1]; }
        _sp_hybrid const float_t& z() const { return (*this)[2]; }

        // Note: we need to delete this constructor or else lambdas
        // with converted omni types <coord> and with <value = float_t>
        // have ambiguous conversion overloads
        point_t(const float_t&) = delete;
    };

    template <typename float_t> std::ostream& operator << (std::ostream& os, const point_t<float_t>& point)
    {
        os << "point: {x: " << point.x() << ", y: " << point.y() << ", z: " << point.z() << "}";
        return os;
    }
}