#include "spade.h"

typedef double real_t;

//we will use a 3x3 variable as an example
template <const int i0, const int i1> using static_dim_t = spade::grid::static_dim_t<i0, i1>;
template <typename integ_t> using dynamic_dim_t = spade::grid::dynamic_dim_t<integ_t>;
using var_map_t   = spade::grid::regular_map_t<static_dim_t<0,3>, static_dim_t<0,3>>;
using block_map_t = spade::grid::regular_map_t<dynamic_dim_t<int>>;
using ijk_map_t   = spade::grid::regular_map_t<dynamic_dim_t<int>, dynamic_dim_t<int>, dynamic_dim_t<int>>;
template <const std::size_t sz> using iv = spade::ctrs::array<int, sz>;
using uni_map_t   = spade::grid::regular_map_t<static_dim_t<0,3>>;
int main(int argc, char** argv)
{
    std::vector<real_t> vec;
    
    var_map_t   vmap;
    block_map_t bmap;
    ijk_map_t   imap;
    
    std::get<0>(bmap.dims) = spade::grid::dynamic_dim_t<int>(0, 8);
    
    std::get<0>(imap.dims) = spade::grid::dynamic_dim_t<int>(-2, 6);
    std::get<1>(imap.dims) = spade::grid::dynamic_dim_t<int>(-2, 6);
    std::get<2>(imap.dims) = spade::grid::dynamic_dim_t<int>(-2, 6);
    imap.compute_coeffs();
    
    spade::ctrs::array<int, 2> i(2,1);
    // print(vmap.offset(i));

    spade::ctrs::array<int, 1> b(2);
    // print(bmap.offset(b));
    
    spade::ctrs::array<int, 3> j(0, 0, 0);
    // print(imap.offset(j));
    
    spade::grid::singleton_map_t nmap;
    spade::grid::composite_map_t cmap(vmap, nmap, imap, nmap, bmap, nmap);
    // uni_map_t a1;
    // spade::grid::composite_map_t zmap(vmap, a1);
    // spade::grid::composite_map_t dmap(a1);
    
    iv<2> i0(1, 2);
    iv<3> i1(3, 4, 0);
    iv<1> i2(2);
    
    // fully verbose indexing
    print(cmap.offset(i0[0], i0[1], i1[0], i1[1], i1[2], i2[0]));
    
    // print(dmap.offset(1));
    // partially verbose indexing (all valid)
    print(cmap.offset(i0,           i1[0], i1[1], i1[2], i2[0]));
    print(cmap.offset(i0[0], i0[1], i1,                  i2[0]));
    print(cmap.offset(i0,           i1[0], i1[1], i1[2], i2));
    print(cmap.offset(i0,           i1,                  i2[0]));
    print(cmap.offset(i0,           i1[0], i1[1], i1[2], i2));
    print(cmap.offset(i0[0], i0[1], i1,                  i2)); 
    
    // compact indexing
    print(cmap.offset(i0, i1, i2));
    // print(cmap.size());
    
    //tiling later
    // using tiled_map_t = spade::grid::tiled_map_t<dynamic_map_t<int>, dynamic_map_t<int>, dynamic_map_t<int>>
    
    // using total_map_t = spade::grid::composite_map_t<var_map_t, ijk_map_t, block_map_t>;
    
    return 0;
}
