#pragma once

#include <iostream>

namespace spade
{
    template <typename dtype, const size_t ar_size> struct bound_box_t
    {
        dtype bnds[2*ar_size];
        
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
        dtype volume(void) const
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