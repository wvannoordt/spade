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
    using cell_t     = spade::grid::cell_idx_t;

    auto lam = [](const point_type& x, const cell_t& i)
    {
        spade::fluid_state::prim_t<real_t> output;
        output.p() = 1 *x[0] + 2 *x[1] + 3 *x[2] + 10*i.lb();
        output.T() = 4 *x[0] + 5 *x[1] + 6 *x[2];
        output.u() = 7 *x[0] + 8 *x[1] + 9 *x[2];
        output.v() = 10*x[0] + 11*x[1] + 12*x[2];
        output.w() = 13*x[0] + 14*x[1] + 15*x[2];
        return output;
    };

    spade::algs::fill_array(prim, lam);
    
    using array_t = decltype(prim);
    {
        using s0_t = spade::omni::prefab::mono_t<spade::grid::cell_centered, spade::omni::info::coord>;
        using s1_t = spade::omni::prefab::mono_t<spade::grid::agno_centered, spade::omni::info::value>;
        using u_t  = spade::omni::stencil_union<s0_t, s1_t>;

        using du_t = spade::omni::stencil_data_t<u_t, array_t>;
        spade::grid::cell_idx_t ic(2,2,2,0);

        du_t du;
        spade::omni::retrieve(prim, ic, du);

        auto du_0 = spade::omni::interpret_stencil<s0_t>(du);
        print(spade::omni::access<spade::omni::info::value>(du_0.cell(0_c)));
        print(prim.get_elem(ic));

        auto du_1 = spade::omni::interpret_stencil<s1_t>(du);
        print(spade::omni::access<spade::omni::info::coord>(du_1.root()));
        print(grid.get_coords(ic));
    }

    {
        using s0_t = spade::omni::prefab::mono_t<spade::grid::agno_centered, spade::omni::info::gradient>;
        using s1_t = spade::omni::prefab::lr_t  <spade::omni::info::value>;
        using u_t  = spade::omni::stencil_union<s0_t, s1_t>;

        using du_t = spade::omni::stencil_data_t<u_t, array_t>;
        spade::grid::face_idx_t ic(0,2,2,2,0);

        du_t du;
        spade::omni::retrieve(prim, ic, du);

        auto du_0 = spade::omni::interpret_stencil<s0_t>(du);
        print(spade::omni::access<spade::omni::info::gradient>(du_0.root()));

        auto du_1 = spade::omni::interpret_stencil<s1_t>(du);
        print(spade::omni::access<spade::omni::info::value>(du_1.cell(0_c)), prim.get_elem(spade::grid::face_to_cell(ic, 0)));
        print(spade::omni::access<spade::omni::info::value>(du_1.cell(1_c)), prim.get_elem(spade::grid::face_to_cell(ic, 1)));
    }

    return 0;
}
