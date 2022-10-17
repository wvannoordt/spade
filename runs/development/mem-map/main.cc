#include <chrono>
#include "spade.h"

typedef double real_t;

//we will use a 3x3 variable as an example
template <const int i0, const int i1> using static_dim_t = spade::grid::static_dim_t<i0, i1>;
using var_map_t   = spade::grid::regular_map_t<static_dim_t<0,3>, static_dim_t<0,3>>;

int main(int argc, char** argv)
{
    std::vector<real_t> vec;
    
    var_map_t vmap;
    // print(vmap.offset());
    
    //tiling later
    // using tiled_map_t = spade::grid::tiled_map_t<dynamic_map_t<int>, dynamic_map_t<int>, dynamic_map_t<int>>
    
    // using ijk_map_t   = spade::grid::composite_map_t<dynamic_map_t<int>, dynamic_map_t<int>, dynamic_map_t<int>>;
    // using block_map_t = spade::grid::dynamic_map_t<int>;
    // using total_map_t = spade::grid::composite_map_t<var_map_t, ijk_map_t, block_map_t>;
    
    return 0;
}
