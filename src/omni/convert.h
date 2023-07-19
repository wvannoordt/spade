#pragma once

#include <concepts>

#include "omni/access.h"

namespace spade::omni
{
    template <typename T> concept has_omni_interface = requires(T t)
    {
        //this may not be the best way to specify this
        typename T::omni_type;
    };

    template <typename list_t, typename kernel_t, typename info_data_t, typename... extracts_t>
    requires(sizeof...(extracts_t) == list_t::num_infos())
    static auto invoke_call(const list_t&, const kernel_t& kernel, const info_data_t& info_list, const extracts_t&... args)
    {
        return kernel(args...);
    }

    template <typename list_t, typename kernel_t, typename info_data_t, typename... extracts_t>
    static auto invoke_call(const list_t&, const kernel_t& kernel, const info_data_t& info_list, const extracts_t&... args)
    {
        constexpr int idx   = sizeof...(extracts_t);
        using info_type     = list_t;
        using query_type    = info_type::template info_elem<idx>;
        const auto& new_arg = access<query_type>(info_list);
        return invoke_call(list_t(), kernel, info_list, args..., new_arg);
    }

    template <typename kernel_t, typename info_data_t>
    static auto invoke_call(const kernel_t& kernel, const info_data_t& info_list)
    {
        using info_type = typename info_data_t::list_type;
        return invoke_call(info_type(), kernel, info_list);
    }

    template <
        typename array_t,
        typename index_t,
        std::invocable<
            info::coord::array_data_type<array_t, index_t::centering_type()>
            > kernel_t
        >
    auto lamda_info_list(const array_t&, const index_t&, const kernel_t& k)
    {
        return info_list_t<info::coord>();
    }

    template <
        typename array_t,
        typename index_t,
        std::invocable<
            info::coord::array_data_type<array_t, index_t::centering_type()>,
            info::index::array_data_type<array_t, index_t::centering_type()>
            > kernel_t
        >
    auto lamda_info_list(const array_t&, const index_t&, const kernel_t& k)
    {
        return info_list_t<info::coord, info::index>();
    }

    template <
        typename array_t,
        typename index_t,
        std::invocable<
            info::value::array_data_type<array_t, index_t::centering_type()>
            > kernel_t
        >
    auto lamda_info_list(const array_t&, const index_t&, const kernel_t& k)
    {
        return info_list_t<info::value>();
    }

    template <
        typename array_t,
        typename index_t,
        std::invocable<
            info::value::array_data_type<array_t, index_t::centering_type()>,
            info::coord::array_data_type<array_t, index_t::centering_type()>
            > kernel_t
        >
    auto lamda_info_list(const array_t&, const index_t&, const kernel_t& k)
    {
        return info_list_t<info::value, info::coord>();
    }

    template <
        typename array_t,
        typename index_t,
        std::invocable<
            info::value::array_data_type<array_t, index_t::centering_type()>,
            info::coord::array_data_type<array_t, index_t::centering_type()>,
            info::index::array_data_type<array_t, index_t::centering_type()>
            > kernel_t
        >
    auto lamda_info_list(const array_t&, const index_t&, const kernel_t& k)
    {
        return info_list_t<info::value, info::coord, info::index>();
    }

    //converts a lambda to an omni object
    template <typename kernel_t, typename array_t, typename index_t>
    struct omni_lambda_wrapper_t
    {
        const kernel_t& kernel;
        using info_list_type = decltype(lamda_info_list(array_t(), index_t(), kernel));
        using omni_type = stencil_t<index_t::centering_type(), elem_t<offset_t<0, 0, 0>, info_list_type>>;
        omni_lambda_wrapper_t(const kernel_t& kernel_in, const array_t&, const index_t&) : kernel{kernel_in}{}
        _sp_hybrid auto operator() (const auto& input_data) const
        {
            const auto elem = input_data.template seek_element<index_t::centering_type()>(0_c);
            return invoke_call(kernel, elem);
        }
    };

    template <const grid::array_centering ctr, typename kernel_t, typename array_t> static auto to_omni(const kernel_t& func, const array_t& array)
    {
        using index_type = grid::get_index_type<ctr>::array_type;
        index_type i = 0;
        if constexpr (omni::has_omni_interface<kernel_t>) {return func;}
        else { return omni::omni_lambda_wrapper_t(func, array, i); }
    }
}