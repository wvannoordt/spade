#pragma once

#include "core/ctrs.h"
#include "omni/omni.h"

namespace spade::convective
{
    template <typename gas_t>
    struct rusanov_t
    {
        const gas_t gas;
        using g_info_type   = typename gas_t::info_type;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = omni::info_union<own_info_type, g_info_type>;
        using float_t       = typename gas_t::value_type;
        using flux_t        = fluid_state::flux_t<float_t>;
        using flux_type     = flux_t;
        rusanov_t(const gas_t& gas_in) : gas{gas_in} {}
        _sp_hybrid ctrs::array<flux_t, 2> operator() (const auto& info) const
        {
            ctrs::array<flux_t, 2> out;
            auto& f_u = out[0];
            auto& f_d = out[1];
            const auto& q              = omni::access<omni::info::value >(info);
            const auto& nv             = omni::access<omni::info::metric>(info);
            const float_t u_n          = nv[0]*q.u() + nv[1]*q.v() + nv[2]*q.w();
            const float_t u_norm       = sqrt(q.u()*q.u() + q.v()*q.v() + q.w()*q.w());
            const float_t rho          = q.p()/(gas.get_R(info)*q.T());
            const float_t pdiff        = q.p()/rho;
            const float_t engy         = float_t(0.5)*(q.u()*q.u()+q.v()*q.v()+q.w()*q.w()) + q.p()/(rho*(gas.get_gamma(info)-float_t(1.0)));
            f_u.continuity()           = float_t(0.5)*rho*u_n;
            f_u.energy()               = float_t(0.5)*rho*u_n*(engy + pdiff);
            f_u.x_momentum()           = float_t(0.5)*rho*q.u()*u_n + float_t(0.5)*q.p()*nv[0];
            f_u.y_momentum()           = float_t(0.5)*rho*q.v()*u_n + float_t(0.5)*q.p()*nv[1];
            f_u.z_momentum()           = float_t(0.5)*rho*q.w()*u_n + float_t(0.5)*q.p()*nv[2];
            f_d = f_u;
            const float_t sigma = float_t(0.5)*(u_norm + sqrt(gas.get_gamma(info)*gas.get_R(info)*q.T()));
            
            f_d.continuity() -= sigma*rho;
            f_d.energy()     -= sigma*rho*engy;
            f_d.x_momentum() -= sigma*rho*q.u();
            f_d.y_momentum() -= sigma*rho*q.v();
            f_d.z_momentum() -= sigma*rho*q.w();
            
            f_u.continuity() += sigma*rho;
            f_u.energy()     += sigma*rho*engy;
            f_u.x_momentum() += sigma*rho*q.u();
            f_u.y_momentum() += sigma*rho*q.v();
            f_u.z_momentum() += sigma*rho*q.w();
            
            return out;
        }
    };

	template <typename ptype>
    struct rusanov_chem_t
    {
        const fluid_state::multicomponent_gas_t<ptype> gas;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = own_info_type;
        using float_t       = ptype;
        using flux_t        = fluid_state::flux_chem_t<float_t>;
		using flux_type     = flux_t;
        rusanov_chem_t(const fluid_state::multicomponent_gas_t<ptype>& gas_in) : gas{gas_in} {}
        _sp_hybrid ctrs::array<flux_t, 2> operator() (const auto& info) const
        {
            ctrs::array<flux_t, 2> out;
            auto& f_u = out[0];
            auto& f_d = out[1];
            const auto& q              = omni::access<omni::info::value >(info);
            const auto& nv             = omni::access<omni::info::metric>(info);

			// Contravariant velocity
            const float_t u_n          = nv[0]*q.u() + nv[1]*q.v() + nv[2]*q.w();

			// Velocity magnitude
            const float_t u_norm       = sqrt(q.u()*q.u() + q.v()*q.v() + q.w()*q.w());

			// Species/Mixture density
			const spade::ctrs::array<float_t, q.ns> rhos = fluid_state::get_rhos(q, gas);
            float_t rho                = 0.0;
			for (int s = 0; s<q.ns; ++s) rho += rhos[s];

			// Vibrational energy
			const float_t Ev           = fluid_state::get_Ev(q, gas);

			// Total Energy
			float_t Etot               = float_t(0.5) * rho * u_norm * u_norm + Ev;
			for (int s = 0; s<q.ns; ++s) Etot += rhos[s] * (gas.get_cvtr(s) * q.T() + gas.hf_s[s]); // Vib energy already added

			// Physical flux
			for (int s = 0; s<q.ns; ++s) f_u.continuity(s) = float_t(0.5) * rhos[s] * u_n;
            f_u.x_momentum()           = float_t(0.5)*rho*q.u()*u_n + float_t(0.5)*q.p()*nv[0];
            f_u.y_momentum()           = float_t(0.5)*rho*q.v()*u_n + float_t(0.5)*q.p()*nv[1];
            f_u.z_momentum()           = float_t(0.5)*rho*q.w()*u_n + float_t(0.5)*q.p()*nv[2];
			f_u.energy()               = float_t(0.5)*u_n*(Etot + q.p());
			f_u.energyVib()            = float_t(0.5)*u_n*Ev;
            f_d = f_u;
            const float_t sigma = float_t(0.5)*(u_norm + fluid_state::get_sos(q, gas));
			
			// Add spectral radius X conservative state vector (upwind)
			for (int s = 0; s<q.ns; ++s) f_u.continuity(s) += sigma * rhos[s];
            f_u.x_momentum() += sigma * rho * q.u();
            f_u.y_momentum() += sigma * rho * q.v();
            f_u.z_momentum() += sigma * rho * q.w();
			f_u.energy()     += sigma * Etot;
			f_u.energyVib()  += sigma * Ev;

			// Add spectral radius X conservative state vector (downwind)
			for (int s = 0; s<q.ns; ++s) f_d.continuity(s) -= sigma * rhos[s];
            f_d.x_momentum() -= sigma * rho * q.u();
            f_d.y_momentum() -= sigma * rho * q.v();
            f_d.z_momentum() -= sigma * rho * q.w();
			f_d.energy()     -= sigma * Etot;
			f_d.energyVib()  -= sigma * Ev;
			
            // Memory address of "out" mapped to f_u and f_d
            return out;
        }
    };
    
    template <typename gas_t>
    struct rusanov_fds_t
    {
        const gas_t gas;
        using g_info_type   = typename gas_t::info_type;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = omni::info_union<own_info_type, g_info_type>;
        using float_t       = typename gas_t::value_type;
        using flux_t        = fluid_state::flux_t<float_t>;
        using flux_type     = flux_t;
        rusanov_fds_t(const gas_t& gas_in) : gas{gas_in} {}
        _sp_hybrid flux_t operator() (const auto& infoF, const auto& qL,const auto& qR) const
        {
            flux_t out{};
            auto& flx = out;
            const auto& qAve           = omni::access<omni::info::value >(infoF);
            const auto& nv             = omni::access<omni::info::metric>(infoF);
            const float_t u_n          = nv[0]*qAve.u() + nv[1]*qAve.v() + nv[2]*qAve.w();
            const float_t rhoL         = qL.p()/(gas.get_R(infoF)*qL.T());
            const float_t engyL        = float_t(0.5)*(qL.u()*qL.u()+qL.v()*qL.v()+qL.w()*qL.w()) + qL.p()/(rhoL*(gas.get_gamma(infoF)-float_t(1.0)));
            const float_t rhoR         = qR.p()/(gas.get_R(infoF)*qR.T());
            const float_t engyR        = float_t(0.5)*(qR.u()*qR.u()+qR.v()*qR.v()+qR.w()*qR.w()) + qR.p()/(rhoR*(gas.get_gamma(infoF)-float_t(1.0)));


            const float_t sigma = (u_n + sqrt(gas.get_gamma(infoF)*gas.get_R(infoF)*qAve.T()));

            flx.continuity() = sigma*rhoL;
            flx.energy()     = sigma*rhoL*engyL;
            flx.x_momentum() = sigma*rhoL*qL.u();
            flx.y_momentum() = sigma*rhoL*qL.v();
            flx.z_momentum() = sigma*rhoL*qL.w();

            flx.continuity() -= sigma*rhoR;
            flx.energy()     -= sigma*rhoR*engyR;
            flx.x_momentum() -= sigma*rhoR*qR.u();
            flx.y_momentum() -= sigma*rhoR*qR.v();
            flx.z_momentum() -= sigma*rhoR*qR.w();

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
        using flux_type     = flux_t;
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
