#pragma once

#include <iostream>
#include "core/range.h"

namespace spade
{
    template <typename dtype, const std::size_t ar_size> struct bound_box_t
    {
        dtype bnds[2*ar_size];
        
        dtype* begin() noexcept {return &bnds[0];}
        dtype* end()   noexcept {return &bnds[0]+2*ar_size;}
        
        bound_box_t(){}
        bound_box_t(const dtype& v)
        {
            for (int i = 0; i < 2*ar_size; ++i) bnds[i] = v;
        }
        const dtype& min(size_t idx) const {return bnds[2*idx+0];}
        const dtype& max(size_t idx) const {return bnds[2*idx+1];}
        
        dtype& min(size_t idx) {return bnds[2*idx+0];}
        dtype& max(size_t idx) {return bnds[2*idx+1];}
        
        dtype size(size_t idx) const {return max(idx)-min(idx);}
        dtype volume() const
        {
            dtype output = 1;
            for (std::size_t i = 0; i < ar_size; i++)
            {
                output *= this->size(i);
            }
            return output;
        }
        
        dtype& operator() (const std::size_t& dim, const std::size_t& min_max) {return bnds[2*dim+min_max];}
        const dtype& operator() (const std::size_t& dim, const std::size_t& min_max) const {return bnds[2*dim+min_max];}

        bound_box_t operator ! () const
        {
            static_assert(std::same_as<dtype, bool>, "cannot apply unary operator ! to non-boolean bound box");
            bound_box_t output = (*this);
            for (auto& b: output) b = !b;
            return output;
        }
    };
    
    template <typename dtype, const size_t ar_size> static std::ostream & operator<<(std::ostream & os, const bound_box_t<dtype, ar_size> & pos)
    {
        os << "{";
        for (auto i: range(0, ar_size))
        {
            os << "(" << pos.min(i) << ", " << pos.max(i) << ")";
            if (i < ar_size-1) os << ", ";
        }
        os << "}";
        return os;
    }
}