#pragma once

#include <tuple>
#include <type_traits>
#include <concepts>
#include <typeinfo>

#include "core/grid.h"
#include "core/coord_system.h"

namespace spade::fetch
{
    template <typename data_t, typename derived_t> struct cf_info_base_t
    {
        data_t data;
    };
    
    template <typename data_t, typename derived_t>
    static std::ostream & operator<<(std::ostream & os, const cf_info_base_t<data_t, derived_t>& ftch)
    {
        os << "Info Type-ID: " << typeid(derived_t()).name() << "\n";
        os << "Data Type-ID: " << typeid(data_t()).name() << "\n";
        os << "Data Value:   " << ftch.data << "\n";
        return os;
    }
    
    template <typename data_t> struct cell_state       : public cf_info_base_t<data_t, cell_state      <data_t>>{};
    template <typename data_t> struct cell_normal      : public cf_info_base_t<data_t, cell_normal     <data_t>>{};
    template <typename data_t> struct cell_index       : public cf_info_base_t<data_t, cell_index      <data_t>>{};
    template <typename data_t> struct cell_coord       : public cf_info_base_t<data_t, cell_coord      <data_t>>{};
    template <typename data_t> struct face_state       : public cf_info_base_t<data_t, face_state      <data_t>>{};
    template <typename data_t> struct face_state_grad  : public cf_info_base_t<data_t, face_state_grad <data_t>>{};
    template <typename data_t> struct face_normal      : public cf_info_base_t<data_t, face_normal     <data_t>>{};
}