#include <spade.h>

#include "k1.h"

template <typename thing_t> void print_type(const thing_t& t)
{
    std::string pf(__PRETTY_FUNCTION__);
    std::size_t start = std::string("void print_type(const thing_t&) [with thing_t = ").length();
    std::size_t end = pf.length()-1;
    print(pf.substr(start, end-start));
}

void print_cen(const spade::grid::array_centering& ac)
{
    switch (ac)
    {
        case spade::grid::edge_centered: {print("edge"); break;}
        case spade::grid::node_centered: {print("node"); break;}
        case spade::grid::face_centered: {print("face"); break;}
        case spade::grid::cell_centered: {print("cell"); break;}
    }
}

int main(int argc, char** argv)
{
    using real_t = double;
    spade::parallel::mpi_t group(&argc, &argv);
    spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    spade::ctrs::array<int, 3> cells_in_block(12, 12, 12);
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
    spade::fluid_state::prim_t<real_t> fill = 0.0;
    fill.p() = 101.315;
    fill.T() = 302.13;
    fill.u() = 65.4;
    fill.v() = 0.005;
    fill.w() = 0.134;
    spade::grid::grid_array prim (grid, fill);
    
    spade::algs::fill_array(prim, [&](){return fill;});
    
    local::kernel0_t k0;

    using stencil_t = local::kernel0_t::stencil_type;
    using array_t   = decltype(prim);
    using data0_t    = spade::omni::stencil_data_t<stencil_t, array_t>;
    
    using data1_t    = spade::omni::stencil_data_t<local::kernel1_t::stencil_type, array_t>;
    
    data0_t p;
    print_type(p.data.data);
    print_type(p.data.next.data);
    print("--------------------");
    print_cen(data0_t::centering_at_i0_elem);
    print_cen(data0_t::next_type::centering_at_i0_elem);
    print_cen(data0_t::next_type::next_type::centering_at_i0_elem);
    print("========================================");
    print_cen(data1_t::centering_at_i0_elem);
    print_cen(data1_t::next_type::centering_at_i0_elem);
    print_cen(data1_t::next_type::next_type::centering_at_i0_elem);
    print_cen(data1_t::next_type::next_type::next_type::centering_at_i0_elem);
    print_cen(data1_t::next_type::next_type::next_type::next_type::centering_at_i0_elem);
    print_cen(data1_t::next_type::next_type::next_type::next_type::next_type::centering_at_i0_elem);
    print_cen(data1_t::next_type::next_type::next_type::next_type::next_type::next_type::centering_at_i0_elem);
    print_cen(data1_t::next_type::next_type::next_type::next_type::next_type::next_type::next_type::centering_at_i0_elem);
    
    data1_t w;
    print("========================================");
    print_cen(spade::omni::detail::ctable_inst_t::elem<0>);
    print_cen(spade::omni::detail::ctable_inst_t::elem<1>);
    print_cen(spade::omni::detail::ctable_inst_t::elem<2>);
    print_cen(spade::omni::detail::ctable_inst_t::elem<3>);
    print_cen(spade::omni::detail::ctable_inst_t::elem<4>);
    print_cen(spade::omni::detail::ctable_inst_t::elem<5>);
    
    using dt2_t = spade::omni::stencil_data_t<local::kernel2_t::stencil_type, array_t>;
    print("========================================");
    print_cen(dt2_t::centering_at_i0_elem);
    print_cen(dt2_t::next_type::centering_at_i0_elem);
    print_cen(dt2_t::next_type::next_type::centering_at_i0_elem);
    print_cen(dt2_t::next_type::next_type::next_type::centering_at_i0_elem);
    print_cen(dt2_t::next_type::next_type::next_type::next_type::centering_at_i0_elem);
    
    print();
    print("Elem counts");
    print("===== s0");
    print("f =", local::kernel0_t::stencil_type::num_face());
    print("c =", local::kernel0_t::stencil_type::num_cell());
    print("n =", local::kernel0_t::stencil_type::num_node());
    
    print("===== s1");
    print("f =", local::kernel1_t::stencil_type::num_face());
    print("c =", local::kernel1_t::stencil_type::num_cell());
    print("n =", local::kernel1_t::stencil_type::num_node());
    
    print("===== s2");
    print("f =", local::kernel2_t::stencil_type::num_face());
    print("c =", local::kernel2_t::stencil_type::num_cell());
    print("n =", local::kernel2_t::stencil_type::num_node());
    
    dt2_t z;
    auto& kk = z.cell(0_c);
    print_type(kk);
    print_type(spade::omni::access<spade::omni::info::value>(kk));
    print_type(spade::omni::access<spade::omni::info::metric>(kk));
    
    auto& ee = spade::omni::access<spade::omni::info::metric>(kk);
    ee.fill(33.03);
    
    // print(spade::omni::access<spade::omni::info::value>(kk));
    print(spade::omni::access<spade::omni::info::metric>(kk));
    
    const spade::grid::cell_idx_t ii(3, 3, 3, 0);
    print(prim.get_elem(ii));
    print(ii);
    
    // ++p;
    
    return 0;
}
