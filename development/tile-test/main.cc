#include <spade.h>

namespace spade::grid
{
    
}

int main(int argc, char** argv)
{
    spade::mem_map::static_dim_t<0, 1> m_vr;
    spade::mem_map::dynamic_dim_t m_x(-2, 22);
    spade::mem_map::dynamic_dim_t m_y(-2, 22);
    spade::mem_map::dynamic_dim_t m_z(-2, 22);
    spade::mem_map::dynamic_dim_t m_lb(0, 12);
    spade::mem_map::tile_uni_dim_t<2, spade::mem_map::dynamic_dim_t<int>> uni_x(m_x);
    spade::mem_map::tile_uni_dim_t<2, spade::mem_map::dynamic_dim_t<int>> uni_y(m_y);
    spade::mem_map::tile_uni_dim_t<2, spade::mem_map::dynamic_dim_t<int>> uni_z(m_z);
    spade::mem_map::tiled_dim_t m_ijk(uni_x, uni_y, uni_z);
    // spade::grid::recti_view_t rct(m_vr, m_ijk, m_lb);
    spade::mem_map::recti_view_t rct(m_ijk, m_lb);
    spade::mem_map::mem_map_t map(rct);
    
    using tile_type = decltype(m_ijk);
    // print(rct.get_extent_array());
    // print(map.i_coeff);
    // print(tile_type::tile_size<0>, tile_type::tile_size<1>, tile_type::tile_size<2>);
    // print(spade::grid::get_dim_from_map<0>(map.map).num_coeffs());
    // print(spade::grid::get_dim_from_map<1>(map.map).num_coeffs());
    // print(spade::grid::get_dim_from_map<2>(map.map).num_coeffs());
    // print(spade::grid::get_dim_from_map<3>(map.map).num_coeffs());
    // print(spade::grid::get_dim_from_map<4>(map.map).num_coeffs());
    spade::grid::cell_idx_t ii(-2,-2,-2,0);
    auto ofst = map.compute_offset(ii);
    print(ofst);
    return 0;
}
