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
    
    using val_t = spade::omni::info::value;
    using mtr_t = spade::omni::info::metric;
    using grd_t = spade::omni::info::gradient;
    using nrm_t = spade::omni::info::normal;

    //stencil 1:
    //
    //
    //             c0    c1
    //     |  +  |  +  |  +  |  +  |
    //
    //
    //
    using o1_t = spade::omni::stencil_t<
        spade::grid::face_centered,
        spade::omni::elem_t<
            spade::omni::offset_t<-1,0,0>,
            spade::omni::info_list_t<val_t,mtr_t,nrm_t>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t< 1,0,0>,
            spade::omni::info_list_t<val_t,mtr_t,nrm_t>
        >
    >;

    //stencil 2:
    //
    //                f0
    //       c0    c1    c2    c3
    //     |  +  |  +  |  +  |  +  |
    //
    //
    //
    using o2_t = spade::omni::stencil_t<
        spade::grid::face_centered,
        spade::omni::elem_t<
            spade::omni::offset_t< 0,0,0>,
            spade::omni::info_list_t<val_t,nrm_t,grd_t>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t<-3,0,0>,
            spade::omni::info_list_t<val_t,mtr_t>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t<-1,0,0>,
            spade::omni::info_list_t<val_t,mtr_t>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t< 1,0,0>,
            spade::omni::info_list_t<val_t,mtr_t>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t< 3,0,0>,
            spade::omni::info_list_t<val_t,mtr_t>
        >
    >;


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
    spade::algs::fill_array(prim, [](const point_type& x, const spade::grid::cell_idx_t& i)
    {
        spade::fluid_state::prim_t<real_t> output;
        output.p() = real_t(i.i());
        output.T() = real_t(i.j());
        output.u() = real_t(i.k());
        output.v() = real_t(i.lb());
        output.w() = 1.0;
        return output;
    });

    //The two stencils disagree on the definition of c0, so when we union the stencils, things run amok
    using ou_t = spade::omni::stencil_union<o1_t, o2_t>;

    using array_t = decltype(prim);
    using d1_t    = spade::omni::stencil_data_t<o1_t, array_t>;
    using d2_t    = spade::omni::stencil_data_t<o2_t, array_t>;
    using du_t    = spade::omni::stencil_data_t<ou_t, array_t>;

    // print_type(ou_t());
    print("==========================================");
    print_type(d1_t());
    print("==========================================");
    print_type(d2_t());
    print("==========================================");
    print_type(du_t());
    print("==========================================");

    spade::grid::face_idx_t j;
    j.i()   = 2;
    j.j()   = 2;
    j.k()   = 2;
    j.lb()  = 0;
    j.dir() = 0;
    du_t data;
    spade::omni::retrieve(prim, j, data);

    // After retrieval, we need to reinterpret the data as relative to the original stencil definitions.

    auto dat1 = spade::omni::interpret_stencil<o1_t>(data);
    auto dat2 = spade::omni::interpret_stencil<o2_t>(data);

    print("absolute stencil union:");
    print(spade::omni::access<spade::omni::info::value>(data.cell(0_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(data.cell(0_c)));
    print(spade::omni::access<spade::omni::info::value>(data.cell(1_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(data.cell(1_c)));
    print(spade::omni::access<spade::omni::info::value>(data.cell(2_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(data.cell(2_c)));
    print(spade::omni::access<spade::omni::info::value>(data.cell(3_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(data.cell(3_c)));

    /*
    print("relative to stencil 1:");
    print(spade::omni::access<spade::omni::info::value>(dat1.cell(0_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(dat1.cell(0_c)));
    print(spade::omni::access<spade::omni::info::value>(dat1.cell(1_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(dat1.cell(1_c)));

    print("relative to stencil 2:");
    print(spade::omni::access<spade::omni::info::value>(dat2.cell(2_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(dat2.cell(2_c)));
    print(spade::omni::access<spade::omni::info::value>(dat2.cell(3_c)), "@ addr", &spade::omni::access<spade::omni::info::value>(dat2.cell(3_c)));
    */
    
    return 0;
}
