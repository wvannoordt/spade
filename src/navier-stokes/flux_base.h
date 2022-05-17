#pragma once

#include <tuple>

namespace cvdf::flux_base
{
    template <typename data_t> struct cell_state
    {
        data_t data;
    };
    
    template <typename data_t> struct cell_normal
    {
        data_t data;
    };
    
    template <typename data_t> struct cell_index
    {
        data_t data;
    };
    
    template <typename data_t> struct cell_coord
    {
        data_t data;
    };
    
    template <typename... infos_t> struct cell_info
    {
        std::tuple<infos_t...> data;
    };
    
    template <typename... infos_t> struct face_info
    {
        std::tuple<infos_t...> data;
    };
    
    template <typename... infos_t> struct left_right
    {
        cell_info<infos_t...> left, right;
    };
    
    template <typename cell_stencil_t, typename face_info_t> struct flux_input_t
    {
        cell_stencil_t cell_data;
        face_info_t    face_data;
    };
    
    template <typename derived_t, typename output_t, typename input_t> struct flux_base_t
    {
        virtual output_t calc_flux(const input_t& flux_input) const {};
    };
}