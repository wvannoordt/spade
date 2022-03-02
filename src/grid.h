#pragma once

#include <concepts>
#include <vector>

#include "output.h"
#include "ctrs.h"
#include "typedef.h"
#include "bounding_box.h"

namespace cvdf::grid
{
    template <class T, typename dtype> concept multiblock_grid = requires(T t, size_t i, size_t j, size_t k, size_t lb)
    {
        // todo: write this
        { t.node_coords(i, j, k, lb) } -> ctrs::vec_nd<3, dtype>;
    };
    
    template <typename dtype> class cartesian_grid_t
    {
        public:
            
            ctrs::array<dtype, 3> node_coords(size_t i, size_t j, size_t k, size_t lb) const
            {
                return ctrs::array<dtype, 3>
            }
            
        private:
            ctrs::array<size_t, dtype> dx;
            ctrs::array<size_t, cvdf_dim> num_blocks;
            ctrs::array<size_t, cvdf_dim> cells_in_block;
            std::vector<bound_box_t<dtype, 3>> block_boxes;
    };
}