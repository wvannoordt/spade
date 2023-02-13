#pragma once

#include <spade.h>

namespace local
{
    struct kernel0_t
    {
        using stencil_type = spade::omni::stencil_t
        <
            spade::grid::cell_centered,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-2, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::index,
                    spade::omni::info::gradient
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 0, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::index,
                    spade::omni::info::gradient
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 2, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::index,
                    spade::omni::info::gradient
                >
            >
        >;
    };
    
    struct kernel1_t
    {        
        using info_type = spade::omni::info_list_t
        <
            spade::omni::info::value,
            spade::omni::info::index,
            spade::omni::info::gradient
        >;
        
        template <const int i> using cell_t = spade::omni::elem_t<spade::omni::offset_t<2*i,0,0>,info_type>;
        
        using stencil_type = spade::omni::stencil_t
        <
            spade::grid::cell_centered,
            cell_t<-1>,
            cell_t< 0>,
            cell_t< 1>
        >;
        
        //Note that the type of "input" necessarily cannot just be the type "stencil type", as
        //the stencil contains neither data nor the types that the data is stored in
        double operator() (const auto& input)
        {
            // const auto qu = spade::omni::retrieve<spade::omni::info::value>(input.cell[2_c]);
            // const auto qh = spade::omni::retrieve<spade::omni::info::value>(input.cell[1_c]);
            // const auto ql = spade::omni::retrieve<spade::omni::info::value>(input.cell[0_c]);
            return 0.0;
        }
    };
}