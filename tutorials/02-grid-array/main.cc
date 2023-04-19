#include "spade.h"

int main(int argc, char** argv)
{
    //We will use double precision for this example
    using real_t = double;

    //Create a parallel group
    spade::parallel::mpi_t group(&argc, &argv);
    
    //Create the dimensions of the group
    spade::ctrs::array<int, 2> num_blocks(2, 2);
    spade::ctrs::array<int, 2> cells_in_block(16, 16);
    spade::ctrs::array<int, 2> exchange_cells(2, 2);

    //Define the physical bounding box of the array
    spade::bound_box_t<real_t, 2> bounds;
    bounds.min(0) =  -1.0;
    bounds.max(0) =   1.0;
    bounds.min(1) =  -1.0;
    bounds.max(1) =   1.0;

    //We won't use any fancy coordinate system for this example,
    //just a plain-old uniform grid.
    spade::coords::identity<real_t> coords;  
    spade::grid::cartesian_grid_t grid(
        num_blocks,
        cells_in_block,
        exchange_cells,
        bounds,
        coords,
        group);
    
    //The type of array we create will be deduced from the
    //element format that we give it!
    real_t fill = 0.0;

    //The array 'arr' is created, valid only over
    //the grid we created
    spade::grid::grid_array arr(grid, fill);
    
    //Here, we can define a lambda that
    //acts as a function to fill our array
    auto ini = [&](const spade::coords::point_t<real_t>& x)
    {
        return x[0] + x[1];
    };
    
    //Fill the array with the function we specified
    spade::algs::fill_array(arr, ini);

    //We can output a VTK file for visualization
    //We specify the directory and the base file name
    spade::io::output_vtk("output", "array", arr);
    
    //We can also output a parallel binary format.
    //This format can be read back in again for later
    //use
    spade::io::binary_write("arr.bin", arr);

    //Here's how we read it back in
    spade::io::binary_read("arr.bin", arr);

    //Note that we can also address specific elements of our array:
    arr(2, 2, 0, 0) = 100.0;

    //Another way to do this is with a cell index (preferred):
    spade::grid::cell_idx_t ii(3, 3, 0, 0);
    arr(ii) = -100.0;

    spade::io::output_vtk("output", "array2", arr);
    return 0;
}
