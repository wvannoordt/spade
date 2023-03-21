#include <iostream>
#include <string>
template <typename thing_t> void print_type(const thing_t& t)
{
    //g++ only
    std::string pf(__PRETTY_FUNCTION__);
    std::size_t start = std::string("void print_type(const thing_t&) [with thing_t = ").length();
    std::size_t end = pf.length()-1;
    std::cout << pf.substr(start, end-start) << std::endl;
}
#include <spade.h>

#include "k1.h"

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
    spade::ctrs::array<int, 3> cells_in_block(10, 10, 10);
    spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    spade::bound_box_t<real_t, 3> bounds;
    bounds.min(0) = 0.0;
    bounds.max(0) = 2.0;
    bounds.min(1) = 0.0;
    bounds.max(1) = 2.0;
    bounds.min(2) = 0.0;
    bounds.max(2) = 2.0;

    spade::coords::identity<real_t> coords;
    spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    spade::fluid_state::prim_t<real_t> fill = 0.0;
    spade::grid::grid_array prim (grid, fill);
    
    using point_type = typename decltype(grid)::coord_point_type;
    spade::algs::fill_array(prim, [](const point_type& x)
    {
        spade::fluid_state::prim_t<real_t> output;
        output.T() = 1 *x[0] + 2 *x[1] + 3 *x[2];
        output.p() = 4 *x[0] + 5 *x[1] + 6 *x[2];
        output.v() = 7 *x[0] + 8 *x[1] + 9 *x[2];
        output.u() = 10*x[0] + 11*x[1] + 12*x[2];
        output.w() = 13*x[0] + 14*x[1] + 15*x[2];
        return output;
    });
    
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
    
    print();
    print("===========");
    
    print(spade::static_math::moddiv<-2, 2>::value);
    print("=================");
    print(-7/4, spade::static_math::moddiv<-7, 4>::value, -2);
    print(-6/4, spade::static_math::moddiv<-6, 4>::value, -2);
    print(-5/4, spade::static_math::moddiv<-5, 4>::value, -2);
    print(-4/4, spade::static_math::moddiv<-4, 4>::value, -1);
    print(-3/4, spade::static_math::moddiv<-3, 4>::value, -1);
    print(-2/4, spade::static_math::moddiv<-2, 4>::value, -1);
    print(-1/4, spade::static_math::moddiv<-1, 4>::value, -1);
    print( 0/4, spade::static_math::moddiv< 0, 4>::value,  0);
    print( 1/4, spade::static_math::moddiv< 1, 4>::value,  0);
    print( 2/4, spade::static_math::moddiv< 2, 4>::value,  0);
    print( 3/4, spade::static_math::moddiv< 3, 4>::value,  0);
    print( 4/4, spade::static_math::moddiv< 4, 4>::value,  1);
    print( 5/4, spade::static_math::moddiv< 5, 4>::value,  1);
    print( 6/4, spade::static_math::moddiv< 6, 4>::value,  1);
    print( 7/4, spade::static_math::moddiv< 7, 4>::value,  1);
    print("=================");
    
    auto modtest = [](auto i0v, auto i1v)
    {
        const int i0 = i0v.value;
        const int i1 = i1v.value;
        print();
        print("intrinsic:", i0, "=", i0/i1, "*", i1, "+", i0%i1, "=", i1*(i0/i1)+(i0%i1));
        print("modular:  ", i0, "=",
            spade::static_math::moddiv<i0,i1>::value, "*", i1,
            "+", spade::static_math::mod<i0,i1>::value, "=",
            i1*(spade::static_math::moddiv<i0,i1>::value)+(spade::static_math::mod<i0,i1>::value));
        print();
    };
    modtest(32_c, 3_c);
    modtest(9_c, 5_c);
    modtest(-16_c, 8_c);
    modtest(-5_c, 2_c);
    modtest(-3_c, 9_c);
    
    {
        print("==============");
        const spade::grid::cell_idx_t ii(0, 0, 0, 0);
        using offs0_t = spade::omni::offset_t<-1,  0,  0>;
        using offs1_t = spade::omni::offset_t< 0, -1,  0>;
        using offs2_t = spade::omni::offset_t< 0,  0, -1>;
        const auto jj0 = offs0_t::compute_index(ii);
        const auto jj1 = offs1_t::compute_index(ii);
        const auto jj2 = offs2_t::compute_index(ii);
        print(ii);
        print(jj0, jj1, jj2);
        print(grid.get_coords(ii));
        print(grid.get_coords(jj0));
        print(grid.get_coords(jj1));
        print(grid.get_coords(jj2));
    }
    
    {
        print("==============");
        const spade::grid::face_idx_t ii(2, 0, 0, 0, 0);
        using offs0_t = spade::omni::offset_t<1,  0,  0>;
        const auto jj0 = offs0_t::compute_index(ii);
        print(ii, jj0);
        print(grid.get_coords(ii), grid.get_coords(jj0));
    }


    // spade::io::output_vtk("output", "prim", prim);

    using o4_type      = local::kernel2_t::stencil_type;
    using o4_data_type = spade::omni::stencil_data_t<o4_type, array_t>;
    const spade::grid::face_idx_t iface(0, 0, 0, 0, 0);
    o4_data_type data;
    spade::omni::retrieve(prim, iface, data);
    print(spade::omni::access<spade::omni::info::gradient>(data.face(0_c)));
    return 0;
}
