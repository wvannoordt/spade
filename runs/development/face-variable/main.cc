#include <chrono>
#include "spade.h"

typedef double real_t;

int main(int argc, char** argv)
{
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(4, 4, 4);
    spade::ctrs::array<int, 3> cells_in_block(16, 16, 16);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    
    bounds.min(0) = 0.0;
    bounds.max(0) = 1.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 1.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 1.0;
    
    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    
    real_t fill = 0.0;
    
    spade::grid::cell_array c_alpha(grid, fill);
    spade::grid::face_array f_alpha(grid, fill);
    
    print(c_alpha.mem_map);
    print(f_alpha.mem_map);
    
    auto ini = [&](const spade::ctrs::array<real_t, 3> x) -> real_t
    {
        return std::sin(x[0]) + std::sin(x[1])*std::cos(x[2]);
    };
    
    spade::utils::mtimer_t tmr("fill_c", "fill_f");
    tmr.start("fill_c");
    spade::algs::fill_array(c_alpha, ini);
    tmr.stop("fill_c");
    print(tmr);
    // spade::grid::cell_idx_t ijkc(0, 0, 0, 0);
    // c_alpha(ijkc) = 99.0;
    
    tmr.start("fill_f");
    spade::algs::fill_array(f_alpha, ini);
    tmr.stop("fill_f");
    // spade::grid::face_idx_t ijkf(0, 0, 0, 0, 0);
    // f_alpha(ijkf) = 99.0;
    
    spade::io::output_vtk("output", "c_alpha", c_alpha);
    spade::io::output_vtk("output", "f_alpha", f_alpha);
    return 0;
}
