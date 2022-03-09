#include "cvdf.h"

int main(int argc, char** argv)
{
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> num_blocks(2, 2, 2);
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> cells_in_block(16, 16, 16);
    cvdf::ctrs::array<size_t, cvdf::cvdf_dim> exchange_cells(2, 2, 2);
    cvdf::bound_box_t<double, cvdf::cvdf_dim> bounds;
    bounds.min(0) = 0.0;
    bounds.min(1) = 0.0;
    bounds.min(2) = 0.0;
    bounds.max(0) = 1.0;
    bounds.max(1) = 1.0;
    bounds.max(2) = 1.0;
    
    cvdf::grid::cartesian_grid_t<double> grid(num_blocks, cells_in_block, exchange_cells, bounds);
    
    return 0;
}