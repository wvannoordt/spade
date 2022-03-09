#pragma once

#include <concepts>
#include <vector>

#include "output.h"
#include "ctrs.h"
#include "typedef.h"
#include "bounding_box.h"
#include "range.h"
#include "static_for.h"

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
            
            cartesian_grid_t(
                const ctrs::array<size_t, cvdf_dim>& num_blocks_in,
                const ctrs::array<size_t, cvdf_dim>& cells_in_block_in,
                const ctrs::array<size_t, cvdf_dim>& exchange_cells_in,
                const bound_box_t<dtype,  cvdf_dim>& bounds_in)
            {
                dx = 1.0;
                num_blocks = 1;
                cells_in_block = 1;
                total_blocks = 1;
                bounds.min(2) = 0.0;
                bounds.max(2) = 1.0;
                
                for (std::size_t i = 0; i < cvdf_dim; i++)
                {
                    bounds.max(i) = bounds_in.max(i);
                    bounds.min(i) = bounds_in.min(i);
                    num_blocks[i] = num_blocks_in[i];
                    cells_in_block[i] = cells_in_block_in[i];
                    dx[i] = bounds.size(i) / (num_blocks[i]*cells_in_block[i]);
                    exchange_cells[i] = exchange_cells_in[i];
                    total_blocks *= num_blocks[i];
                }
                block_boxes.resize(total_blocks);
                std::size_t clb = 0;
                for (auto lb: range(0,num_blocks[0])*range(0,num_blocks[1])*range(0,num_blocks[2]))
                {
                    auto& box = block_boxes[clb];
                    ctrs::array<dtype, 3> lower;
                    ctrs::array<dtype, 3> upper;
                    static_for<0,3>([&](auto i)
                    {
                        box.min(i.value) = bounds.min(i.value) + (lb[i.value]+0)*bounds.size(i.value)/num_blocks[i.value];
                        box.max(i.value) = bounds.min(i.value) + (lb[i.value]+1)*bounds.size(i.value)/num_blocks[i.value];
                    });
                    ++clb;
                }
            }
            
            ctrs::array<dtype, 3> node_coords(size_t i, size_t j, size_t k, size_t lb) const
            {
                ctrs::array<dtype, 3> output(0, 0, 0);
                auto& box = block_boxes[lb];
                for (size_t idir = 0; idir < cvdf_dim; idir++)
                {
                    output[idir] = box.min(idir) + i*dx[idir];
                }
                return output;
            }
            
        private:
            ctrs::array<dtype,  3> dx;
            ctrs::array<size_t, 3> num_blocks;
            ctrs::array<size_t, 3> cells_in_block;
            ctrs::array<size_t, 3> exchange_cells;
            bound_box_t<dtype,  3> bounds;
            std::vector<bound_box_t<dtype, 3>> block_boxes;
            size_t total_blocks;
    };
}