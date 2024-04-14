#pragma once

#include "omni/omni.h"

namespace spade::state_sensor
{
    template <typename derived_t, typename ilist_t = omni::info_list_t<>>
    struct sensor_interface_t
    {
        _sp_hybrid derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        _sp_hybrid const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = ilist_t;

        _sp_hybrid auto get_sensor(const auto& input) const
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

        _sp_hybrid ducros_t(const float_t& epsilon_in) : epsilon{epsilon_in} {}

        _sp_hybrid float_t get_sensor(const auto& info) const
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

    template <typename float_t> struct ducros_outramp_t
		: public sensor_interface_t<ducros_outramp_t<float_t>, omni::info_list_t<omni::info::gradient,spade::omni::info::coord>>
    {
        using base_t = sensor_interface_t<ducros_outramp_t<float_t>, omni::info_list_t<omni::info::gradient,spade::omni::info::coord>>;
        using base_t::get_sensor;
        using base_t::info_type;

        float_t epsilon;
		float_t x_out;

        _sp_hybrid ducros_outramp_t(const float_t& x_out_in,const float_t& epsilon_in) : epsilon{epsilon_in},x_out{x_out_in} {}

        _sp_hybrid float_t get_sensor(const auto& info) const
        {
            const auto& grad     = omni::access<omni::info::gradient>(info);
			const auto& x      = spade::omni::access<spade::omni::info::coord>(info);
            auto        theta_sq = grad[0].u() + grad[1].v() + grad[2].w();
            theta_sq *= theta_sq;
            const auto vort0     = grad[1].w() - grad[2].v();
            const auto vort1     = grad[2].u() - grad[0].w();
            const auto vort2     = grad[0].v() - grad[1].u();
            const auto vort_sq   = vort0*vort0 + vort1*vort1 + vort2*vort2;
			auto value = (x[0]>x_out)? (1.0) : (theta_sq / (theta_sq + vort_sq + epsilon));
			return value;
        }
    };

    template <typename float_t> struct ducros_ibm_t
		: public sensor_interface_t<ducros_ibm_t<float_t>, omni::info_list_t<omni::info::gradient>>
    {
        using base_t = sensor_interface_t<ducros_ibm_t<float_t>, omni::info_list_t<omni::info::gradient>>;
        using base_t::get_sensor;
        using base_t::info_type;

        float_t epsilon;
		float_t x_out;
		float_t irreg_alpha;

        _sp_hybrid ducros_ibm_t(const float_t& irreg_alpha_in,const float_t& epsilon_in) : epsilon{epsilon_in},irreg_alpha{irreg_alpha_in} {}

        _sp_hybrid float_t get_sensor(const auto& info) const
        {
            const auto& grad     = omni::access<omni::info::gradient>(info);
			const auto& x      = spade::omni::access<spade::omni::info::coord>(info);
            auto        theta_sq = grad[0].u() + grad[1].v() + grad[2].w();
            theta_sq *= theta_sq;
            const auto vort0     = grad[1].w() - grad[2].v();
            const auto vort1     = grad[2].u() - grad[0].w();
            const auto vort2     = grad[0].v() - grad[1].u();
            const auto vort_sq   = vort0*vort0 + vort1*vort1 + vort2*vort2;
			auto value = max(irreg_alpha,(theta_sq / (theta_sq + vort_sq + epsilon)));
			return value;
        }
    };

    template <typename float_t> struct const_sensor_t
		: public sensor_interface_t<const_sensor_t<float_t>, omni::info_list_t<>>
    {
        using base_t = sensor_interface_t<const_sensor_t<float_t>, omni::info_list_t<>>;
        using base_t::get_sensor;
        using base_t::info_type;

        float_t alpha;

        _sp_hybrid const_sensor_t(const float_t& alpha_in) : alpha{alpha_in} {}

        _sp_hybrid float_t get_sensor(const auto& info) const
        {
            return alpha;
        }
    };

    template <typename float_t> struct position_sensor_t
		: public sensor_interface_t<position_sensor_t<float_t>, omni::info_list_t<spade::omni::info::coord>>
    {
        using base_t = sensor_interface_t<position_sensor_t<float_t>, omni::info_list_t<spade::omni::info::coord>>;
        using base_t::get_sensor;
        using base_t::info_type;

		float_t x_out;
        float_t alpha;

        _sp_hybrid position_sensor_t(const float_t& x_in, const float_t& alpha_in) : x_out{x_in},alpha{alpha_in} {}

        _sp_hybrid float_t get_sensor(const auto& info) const
        {
			const auto& x      = spade::omni::access<spade::omni::info::coord>(info);
            return (x[0]>x_out) ? alpha:float_t(0.);
        }
    };

}
