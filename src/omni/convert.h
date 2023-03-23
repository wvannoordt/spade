#pragma once

#include <concepts>

namespace spade::omni
{
    template <typename T> concept has_omni_interface = requires(T t)
    {
        //this may not be the best way to specify this
        T::omni_type();
    };

    namespace detail
    {
        template <typename index_t, typename kernel_t, typename info_data_t, typename... extracts_t>
        static auto invoke_call(const kernel_t& kernel, const info_data_t& info_list, const extracts_t&... args)
        {
            if constexpr (requires {kernel(args...);})
            {
                return kernel(args...);
            }
            else
            {
                return invoke_call(kernel, info_list, args..., info_list.data);
            }
        }
    }

    template <typename array_t, typename index_t, std::invocable<info::coord::array_data_type<array_t, index_t>> kernel_t> auto type_test(const kernel_t& k) {return info_list_t<info::coord>();}

    //converts a lambda to an omni object
    template <typename kernel_t, typename array_t, typename index_t>
    struct omni_lambda_wrapper_t
    {
        const kernel_t& kernel;
        using omni_type = prefab::mono_t<index_t::centering_type(), decltype(type_test(kernel_t()))>;
        omni_lambda_wrapper_t(const kernel_t& kernel_in, const array_t&, const index_t&) : kernel{kernel_in}{}
        auto operator() (const auto& input_data) const
        {
            const auto elem = input_data.template seek_element<index_t::centering_type()>(0_c);
            return detail::invoke_call(kernel, elem);
        }
    };

    template <typename kernel_t> static auto to_omni(const kernel_t& func)
    {
        if constexpr (omni::has_omni_interface<kernel_t>) {return kernel;}
        else { return omni::omni_lambda_wrapper_t(kernel, array, i); }
    }
}