#pragma once

#include <spade.h>

namespace local
{
    struct kernel0_t
    {
        using stencil_type = spade::omni::stencil_t
        <
            spade::grid::face_centered,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-1, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::metric
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 0, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::gradient,
                    spade::omni::info::metric
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 1, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::metric
                >
            >
        >;
    };
    
    struct kernel1_t
    {
        using stencil_type = spade::omni::stencil_t
        <
            spade::grid::cell_centered,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-1, -1, -1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 1, -1, -1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-1,  1, -1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 1,  1, -1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-1, -1,  1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 1, -1,  1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-1,  1,  1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t< 1,  1,  1>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value
                >
            >
        >;
    };
    
    struct kernel2_t
    {
        using stencil_type = spade::omni::stencil_t
        <
            spade::grid::face_centered,
            spade::omni::elem_t
            <
                spade::omni::offset_t<0, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::gradient,
                    spade::omni::info::metric
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-3, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::metric
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-1, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::metric
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<1, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::metric
                >
            >,
            spade::omni::elem_t
            <
                spade::omni::offset_t<3, 0, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::metric
                >
            >
        >;
    };
    
    struct kernel3_t
    {
        using stencil_type = spade::omni::stencil_t
        <
            spade::grid::cell_centered,
            spade::omni::elem_t
            <
                spade::omni::offset_t<-2, -2, 0>,
                spade::omni::info_list_t
                <
                    spade::omni::info::value,
                    spade::omni::info::gradient,
                    spade::omni::info::metric
                >
            >
        >;
    };
}