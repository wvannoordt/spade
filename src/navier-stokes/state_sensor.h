#pragma once

#include "omni/omni.h"

namespace spade::state_sensor
{
    template <typename derived_t, typename ilist_t = omni::info_list_t<>>
    struct sensor_interface_t
    {
        derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = ilist_t;

        auto get_sensor(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_sensor(args...);}, input);
        }
    };

    template <typename float_t> struct ducros_t
    : public sensor_interface_t<ducros_t<float_t>, omni::info_list_t<omni::info::gradient>>
    {
        using base_t = sensor_interface_t<ducros_t<float_t>, omni::info_list_t<omni::info::gradient>>;
        using base_t::get_sensor;
        using base_t::info_type;

        float_t epsilon;

        ducros_t(const float_t& epsilon_in) : epsilon{epsilon_in} {}

        float_t get_sensor(const auto& info) const
        {
            const auto& grad     = omni::access<omni::info::gradient>(info);
            auto        theta_sq = grad[0].u() + grad[1].v() + grad[2].w();
            theta_sq *= theta_sq;
            const auto vort0     = grad[1].w() - grad[2].v();
            const auto vort1     = grad[2].u() - grad[0].w();
            const auto vort2     = grad[0].v() - grad[1].u();
            const auto vort_sq   = vort0*vort0 + vort1*vort1 + vort2*vort2;
            return theta_sq / (theta_sq + vort_sq + epsilon);
        }
    };
}