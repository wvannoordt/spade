#pragma once

// #ifdef SDF_MAJOR_VERSION
#include "scidf.h"
#include <vector>
#include "core/ctrs.h"
namespace scidf
{
    template <typename data_t, const std::size_t ar_size>
    const scidf::iconversion_t& operator >> (const scidf::iconversion_t& cnv, spade::ctrs::array<data_t, ar_size>& data)
    {
        std::vector<data_t> input;
        cnv >> input;
        if (input.size() != ar_size) throw sdf_exception("incorrect number of elements: expecting " + std::to_string(ar_size) + ", got " + std::to_string(input.size()));
        for (int i = 0; i < ar_size; ++i) data[i] = input[i];
        return cnv;
    }

    template <typename data_t, const std::size_t ar_dim>
    const scidf::iconversion_t& operator >> (const scidf::iconversion_t& cnv, spade::bound_box_t<data_t, ar_dim>& data)
    {
        std::vector<std::vector<data_t>> input;
        // e.g. [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
        cnv >> input;
        if (input.size() != ar_dim) throw sdf_exception("incorrect number of elements: expecting " + std::to_string(ar_dim) + ", got " + std::to_string(input.size()));
        for (int i = 0; i < ar_dim; ++i)
        {
            if (input[i].size() != 2) throw sdf_exception("expecting 2 entries for dimension bounds when parsing bounding box, received " + std::to_string(input[i].size()));
        }

        for (int i = 0; i < ar_dim; ++i)
        {
            data.min(i) = input[i][0];
            data.max(i) = input[i][1];
        }
        return cnv;
    }
}
// #endif