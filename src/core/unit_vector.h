#pragma once

#include "core/cuda_incl.h"

namespace spade::utils
{
    // Probably didn't need a whole header for this
    struct unit_vector_t
    {
        int dir;
        _sp_hybrid int operator[] (int ii) const { return int(ii == dir); }
    };
    
    inline _sp_hybrid unit_vector_t unit(int ii) { return unit_vector_t{ii};}
}