#pragma once

#include "core/grid_index_types.h"
#include "core/utils.h"
#include "core/grid.h"

namespace spade::omni
{
    namespace info
    {
        template <typename T> concept is_omni_info = requires(const T& t) {T::info_tag;};
        template <typename derived_t> struct info_base
        {
            constexpr static bool info_tag = true;
        };
        
        
        struct index    : public info_base<index    >
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
        
        struct normal   : public info_base<normal   > // only supported at locations centered at faces, expressed in index space
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
        
        struct coord    : public info_base<coord    >
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
        
        struct metric   : public info_base<metric   >
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
        
        struct jacobian : public info_base<jacobian >
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
        
        struct value    : public info_base<value    > // interpolation to nodes?
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
        
        struct gradient : public info_base<gradient > // floating-point types only, and need to be able to specify the order-of-accuracy later
        {
            template <typename array_t, const grid::array_centering center>
            using array_data_type = typename grid::get_index_type<center>::array_type;
        };
    }
    template <typename... infos_t>
    requires (info::is_omni_info<infos_t> && ...)
    struct info_list_t
    {
        constexpr static int num_infos() {return sizeof...(infos_t);}
        
        template <const int idx>
        using info_elem = typename utils::get_pack_type<idx, infos_t...>::type;
    };
}