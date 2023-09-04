#pragma once

#include "core/ctrs.h"
#include "omni/omni.h"

namespace spade::convective
{
    template <typename gas_t> struct rusanov_t
    {
        const gas_t gas;
        using g_info_type   = typename gas_t::info_type;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = omni::info_union<own_info_type, g_info_type>;
        using float_t       = typename gas_t::value_type;
        using flux_t        = fluid_state::flux_t<float_t>;
        rusanov_t(const gas_t& gas_in) : gas{gas_in} {}
        _sp_hybrid ctrs::array<flux_t, 2> operator() (const auto& info) const
        {
            ctrs::array<flux_t, 2> out;
            auto& f_u = out[0];
            auto& f_d = out[1];
            const auto& q              = omni::access<omni::info::value >(info);
            const auto& nv             = omni::access<omni::info::metric>(info);
            const float_t u_n          = nv[0]*q.u() + nv[1]*q.v() + nv[2]*q.w();
            const float_t rho          = q.p()/(gas.get_R(info)*q.T());
            const float_t enth         = 0.5*(q.u()*q.u()+q.v()*q.v()+q.w()*q.w()) + q.p()/rho + q.p()/(rho*(gas.get_gamma(info)-1.0));
            f_u.continuity()           = 0.5*rho*u_n;
            f_u.energy()               = 0.5*rho*u_n*enth;
            f_u.x_momentum()           = 0.5*rho*q.u()*u_n + 0.5*q.p()*nv[0];
            f_u.y_momentum()           = 0.5*rho*q.v()*u_n + 0.5*q.p()*nv[1];
            f_u.z_momentum()           = 0.5*rho*q.w()*u_n + 0.5*q.p()*nv[2];
            f_d = f_u;
            const real_t sigma = abs(u_n) + sqrt(gas.get_gamma(info)*gas.get_R(info)*q.T());
            
            //Note: enth is actually the enthalpy as implemented, 
            //should actually be the conserved quantity total energy (pressure diffusion term)
            f_d.continuity() -= sigma*rho;
            f_d.energy()     -= sigma*rho*enth;
            f_d.x_momentum() -= sigma*rho*q.u();
            f_d.y_momentum() -= sigma*rho*q.v();
            f_d.z_momentum() -= sigma*rho*q.w();
            
            f_u.continuity() += sigma*rho;
            f_u.energy()     += sigma*rho*enth;
            f_u.x_momentum() += sigma*rho*q.u();
            f_u.y_momentum() += sigma*rho*q.v();
            f_u.z_momentum() += sigma*rho*q.w();
            
            return out;
        }
    };

    template <typename gas_t> struct phys_flux_t
    {
        const gas_t gas;
        using g_info_type   = typename gas_t::info_type;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = omni::info_union<own_info_type, g_info_type>;
        using float_t       = typename gas_t::value_type;
        using flux_t        = fluid_state::flux_t<float_t>;
        phys_flux_t(const gas_t& gas_in) : gas{gas_in} {}
        _sp_hybrid ctrs::array<flux_t, 2> operator() (const auto& info) const
        {
            ctrs::array<flux_t, 2> out;
            auto& f_u = out[0];
            auto& f_d = out[1];
            const auto& q              = omni::access<omni::info::value >(info);
            const auto& nv             = omni::access<omni::info::metric>(info);
            const float_t u_n          = nv[0]*q.u() + nv[1]*q.v() + nv[2]*q.w();
            const float_t rho          = q.p()/(gas.get_R(info)*q.T());
            const float_t enth         = 0.5*(q.u()*q.u()+q.v()*q.v()+q.w()*q.w()) + q.p()/rho + q.p()/(rho*(gas.get_gamma(info)-1.0));
            f_u.continuity()           = 0.5*rho*u_n;
            f_u.energy()               = 0.5*rho*u_n*enth;
            f_u.x_momentum()           = 0.5*rho*q.u()*u_n + 0.5*q.p()*nv[0];
            f_u.y_momentum()           = 0.5*rho*q.v()*u_n + 0.5*q.p()*nv[1];
            f_u.z_momentum()           = 0.5*rho*q.w()*u_n + 0.5*q.p()*nv[2];
            f_d = f_u;
            return out;
        }
    };
}