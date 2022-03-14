#pragma once

#include <iostream>
#include <fstream>
#include "grid.h"
namespace cvdf::output
{    
    template <class output_stream_t, grid::multiblock_grid grid_output_t, typename... arrays_t>
    void output_grid(output_stream_t& out_str, const grid_output_t& obj, arrays_t... arrays)
    {
        out_str << "hello\n";
    }
}