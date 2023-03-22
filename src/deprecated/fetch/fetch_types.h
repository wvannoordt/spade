#pragma once

#include <tuple>
#include <type_traits>
#include <typeinfo>

#include "core/config.h"
#include "core/grid.h"
#include "core/coord_system.h"

namespace spade::fetch
{
    enum fetch_centering
    {
        fetch_face,
        fetch_cell,
        fetch_node
    };
    
    template <typename cell_stencil_t, typename face_info_t> struct face_fetch_t
    {
        static constexpr fetch_centering centering() {return fetch_face; }
        cell_stencil_t cell_data;
        face_info_t    face_data;
    };

    template <typename cell_stencil_t, typename face_info_t>
    static std::ostream & operator<<(std::ostream & os, const face_fetch_t<cell_stencil_t, face_info_t>& ftch)
    {
        os << "Face fetch:\n >> Cell data:\n" << ftch.cell_data << "\n >> Face data:\n" << ftch.face_data;
        return os;
    }

    template <typename cell_stencil_t> struct cell_fetch_t
    {
        static constexpr fetch_centering centering() {return fetch_cell; }
        cell_stencil_t cell_data;
    };

    template <typename cell_stencil_t>
    static std::ostream & operator<<(std::ostream & os, const cell_fetch_t<cell_stencil_t>& ftch)
    {
        os << "Cell fetch:\n >> Cell data:\n" << ftch.cell_data;
        return os;
    }
    
    template <typename T, const fetch_centering ctr> concept fetch_type = (T::centering()==ctr);
}