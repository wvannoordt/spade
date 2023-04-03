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
    //       c0    c1
    //     |  +  |  +  |
    //
    //
    //
    using o1_t = spade::omni::stencil_t<
        spade::grid::face_centered,
        spade::omni::elem_t<
            spade::omni::offset_t<-1,0,0>,
            spade::omni::info_list_t<val_t,mtr_t>
        >,
        spade::omni::elem_t<
            spade::omni::offset_t< 1,0,0>,
            spade::omni::info_list_t<val_t,mtr_t>
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

    
    return 0;
}
