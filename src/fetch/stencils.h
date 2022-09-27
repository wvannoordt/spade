#pragma once

#include <tuple>
#include <type_traits>
#include <concepts>
#include <typeinfo>

#include "core/grid.h"
#include "core/coord_system.h"

namespace spade::fetch
{
    template <typename cell_info_t> struct left_right
    {
        cell_info_t left, right;
    };
    
    template <typename cell_info_t>
    static std::ostream & operator<<(std::ostream & os, const left_right<cell_info_t>& ftch)
    {
        os << "Left-Right-stencil\nData:\n";
        os << "Left:\n"  << ftch.left  << "\n";
        os << "Right:\n" << ftch.right << "\n";
        return os;
    }
    
    template <const std::size_t stencil_size, typename cell_info_t> struct flux_line
    {
        ctrs::array<cell_info_t, stencil_size> stencil;
    };
    
    template <const std::size_t stencil_size, typename cell_info_t>
    static std::ostream & operator<<(std::ostream & os, const flux_line<stencil_size, cell_info_t>& ftch)
    {
        os << "Line-stencil, size " << stencil_size <<"\nData:\n";
        for (auto i: range(0, stencil_size))
        {
            os << "Stencil location " << i << ":\n" << ftch.stencil[i] << "\n";
        }
        return os;
    }
    
    template <typename cell_info_t> struct cell_mono
    {
        cell_info_t mono_data;
    };
    
    template <typename cell_info_t>
    static std::ostream & operator<<(std::ostream & os, const cell_mono<cell_info_t>& ftch)
    {
        os << "Mono-stencil\nData:\n";
        os << ftch.mono_data << "\n";
        return os;
    }
}