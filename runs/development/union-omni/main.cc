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
    // spade::parallel::mpi_t group(&argc, &argv);
    // spade::ctrs::array<int, 3> num_blocks(2, 2, 2);
    // spade::ctrs::array<int, 3> cells_in_block(10, 10, 10);
    // spade::ctrs::array<int, 3> exchange_cells(2, 2, 2);
    // spade::bound_box_t<real_t, 3> bounds;
    // bounds.min(0) = 0.0;
    // bounds.max(0) = 2.0;
    // bounds.min(1) = 0.0;
    // bounds.max(1) = 2.0;
    // bounds.min(2) = 0.0;
    // bounds.max(2) = 2.0;

    // spade::coords::identity<real_t> coords;
    // spade::grid::cartesian_grid_t grid(num_blocks, cells_in_block, exchange_cells, bounds, coords, group);
    // spade::fluid_state::prim_t<real_t> fill = 0.0;
    // spade::grid::grid_array prim (grid, fill);

    using val_t = spade::omni::info::value;
    using met_t = spade::omni::info::metric;
    using gra_t = spade::omni::info::gradient;
    using nor_t = spade::omni::info::normal;

    using o1_inf_t = spade::omni::info_list_t<val_t, gra_t>;
    using o2_inf_t = spade::omni::info_list_t<gra_t, met_t>;
    using o3_inf_t = spade::omni::info_list_t<val_t, nor_t>;

    using ou_inf_t = spade::omni::info_union<o1_inf_t, o2_inf_t, o3_inf_t>;
    print(o1_inf_t::contains<val_t>, o2_inf_t::contains<val_t>, o3_inf_t::contains<val_t>, "->", ou_inf_t::contains<val_t>);
    print(o1_inf_t::contains<gra_t>, o2_inf_t::contains<gra_t>, o3_inf_t::contains<gra_t>, "->", ou_inf_t::contains<gra_t>);
    print(o1_inf_t::contains<met_t>, o2_inf_t::contains<met_t>, o3_inf_t::contains<met_t>, "->", ou_inf_t::contains<met_t>);
    print(o1_inf_t::contains<nor_t>, o2_inf_t::contains<nor_t>, o3_inf_t::contains<nor_t>, "->", ou_inf_t::contains<nor_t>);


    return 0;
}
