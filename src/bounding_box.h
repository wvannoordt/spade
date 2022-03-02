#pragma once

namespace cvdf
{
    template <typename dtype, const size_t ar_size> struct bound_box_t
    {
        dtype bnds[2*ar_size];
        
        const dtype& min(size_t idx) const {return bnds[2*idx+0];}
        const dtype& max(size_t idx) const {return bnds[2*idx+1];}
        
        dtype& min(size_t idx) {return bnds[2*idx+0];}
        dtype& max(size_t idx) {return bnds[2*idx+1];}
        
        dtype size(size_t idx) {return max(idx)-min(idx);}
    };
}