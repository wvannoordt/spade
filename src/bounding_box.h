#pragma once

#include <iostream>

namespace cvdf
{
    template <typename dtype, const size_t ar_size> struct bound_box_t
    {
        dtype bnds[2*ar_size];
        
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
        
    };
}