#pragma once

#include "core/cuda_incl.h"
#include "omni/omni.h"

namespace spade::fluid_state
{
    template <class T> concept is_state_type = requires(T t, size_t idx)
    {
        t[idx];
        T::size();
        t.name(0);
        typename T::value_type;
    };

    template <typename derived_t, typename ilist_t = omni::info_list_t<>>
    struct gas_interface_t
    {
        derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = ilist_t;

        auto get_R    (const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_R    (args...);}, input);
        }

        auto get_gamma(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_gamma(args...);}, input);
        }
    };

    template <typename dtype> struct ideal_gas_t
    : public gas_interface_t<ideal_gas_t<dtype>>
    {
        using base_t = gas_interface_t<ideal_gas_t<dtype>>;
        using base_t::get_R;
        using base_t::get_gamma;
        using base_t::info_type;

        typedef dtype value_type;
        dtype R, gamma;
        ideal_gas_t(){}
        ideal_gas_t(const dtype& gamma_in, const dtype& R_in) : gamma{gamma_in}, R{R_in} {}
        dtype get_R    () const {return this->R;}
        dtype get_gamma() const {return this->gamma;}
    };
}