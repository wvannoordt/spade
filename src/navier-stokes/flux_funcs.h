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
    
    template <typename gas_t>
    struct rusanov_fds_t
    {
        const gas_t gas;
        using g_info_type   = typename gas_t::info_type;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = omni::info_union<own_info_type, g_info_type>;
        using float_t       = typename gas_t::value_type;
        using flux_t        = fluid_state::flux_t<float_t>;
        rusanov_fds_t(const gas_t& gas_in) : gas{gas_in} {}
        _sp_hybrid flux_t operator() (const auto& infoF, const auto& qL, const auto& qR) const
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

    // HLLC Flux Scheme <JRB | Implemented: 4-12-24 | Validated: TODO>
    template <typename gas_t>
    struct hllc_fds_t
    {
        const gas_t gas;
        using g_info_type   = typename gas_t::info_type;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using info_type     = omni::info_union<own_info_type, g_info_type>;
        using float_t       = typename gas_t::value_type;
        using flux_t        = fluid_state::flux_t<float_t>;
        hllc_fds_t(const gas_t& gas_in) : gas{gas_in} {}
        _sp_hybrid flux_t operator() (const auto& infoF, const auto& qL, const auto& qR) const
        {
            flux_t out{};
            auto& flx = out;
            const auto& nv             = omni::access<omni::info::metric>(infoF);
            const float_t uL_n         = nv[0]*qL.u() + nv[1]*qL.v() + nv[2]*qL.w();
            const float_t uR_n         = nv[0]*qR.u() + nv[1]*qR.v() + nv[2]*qR.w();
            const float_t rhoL         = qL.p()/(gas.get_R(infoF)*qL.T());
            const float_t engyL        = float_t(0.5)*(qL.u()*qL.u()+qL.v()*qL.v()+qL.w()*qL.w()) + qL.p()/(rhoL*(gas.get_gamma(infoF)-float_t(1.0)));
            const float_t rhoR         = qR.p()/(gas.get_R(infoF)*qR.T());
            const float_t engyR        = float_t(0.5)*(qR.u()*qR.u()+qR.v()*qR.v()+qR.w()*qR.w()) + qR.p()/(rhoR*(gas.get_gamma(infoF)-float_t(1.0)));

            const float_t aL            = sqrt(gas.get_gamma(infoF)*gas.get_R(infoF)*qL.T());
            const float_t aR            = sqrt(gas.get_gamma(infoF)*gas.get_R(infoF)*qR.T());
            const float_t z             = (gas.get_gamma(infoF)-1)/(2*gas.get_gamma(infoF));
            const float_t p_star        = std::pow((aL+aR-(gas.get_gamma(infoF)-1)/2*(uR_n-uL_n))/(aL/std::pow(qL.p(), z)+aR/std::pow(qR.p(), z)), 1/z);
            const float_t sL            = uL_n - (p_star <= qL.p())? aL : aL*sqrt(1+(gas.get_gamma(infoF)+1)/(2*gas.get_gamma(infoF))*(p_star/qL.p()-1));
            const float_t sR            = uR_n + (p_star <= qR.p())? aR : aR*sqrt(1+(gas.get_gamma(infoF)+1)/(2*gas.get_gamma(infoF))*(p_star/qR.p()-1));
            const float_t s_star        = (qR.p()-qL.p()+rhoL*uL_n*(sL-uL_n)-rhoR*uR_n*(sR-uR_n))/(rhoL*(sL-uL_n)-rhoR*(sR-uR_n));

            if (sL >= 0 || (sL < 0 && s_star >= 0)) 
            {
                flx.continuity() = rhoL*uL_n;
                flx.x_momentum() = rhoL*qL.u()*uL_n+qL.p()*nv[0];
                flx.y_momentum() = rhoL*qL.v()*uL_n+qL.p()*nv[1];
                flx.z_momentum() = rhoL*qL.w()*uL_n+qL.p()*nv[2];
                flx.energy()     = uL_n*(engyL+qL.p());

                if (sL < 0 && s_star >= 0)
                {
                    const float_t rho_star = rhoL*(sL-uL_n)/(sL-s_star);
                    flx.continuity() += sL*(rho_star - rhoL);
                    flx.x_momentum() += sL*(rho_star*(s_star*nv[0] + qL.u()*(1-nv[0])) - rhoL*qL.u());
                    flx.y_momentum() += sL*(rho_star*(s_star*nv[1] + qL.v()*(1-nv[1])) - rhoL*qL.v());
                    flx.z_momentum() += sL*(rho_star*(s_star*nv[2] + qL.w()*(1-nv[2])) - rhoL*qL.w());
                    flx.energy()     += sL*(rho_star*(engyL/rhoL+(s_star-uL_n)*(s_star+qL.p()/(rhoL*(sL-uL_n)))) - engyL);
                }
            }
            else
            {
                flx.continuity() = rhoR*uR_n;
                flx.x_momentum() = rhoR*qR.u()*uR_n+qR.p()*nv[0];
                flx.y_momentum() = rhoR*qR.v()*uR_n+qR.p()*nv[1];
                flx.z_momentum() = rhoR*qR.w()*uR_n+qR.p()*nv[2];
                flx.energy()     = uR_n*(engyR+qR.p());

                if (sR > 0 && s_star <= 0)
                {
                    const float_t rho_star = rhoR*(sR-uR_n)/(sR-s_star);
                    flx.continuity() += sR*(rho_star - rhoR);
                    flx.x_momentum() += sR*(rho_star*(s_star*nv[0] + qR.u()*(1-nv[0])) - rhoR*qR.u());
                    flx.y_momentum() += sR*(rho_star*(s_star*nv[1] + qR.v()*(1-nv[1])) - rhoR*qR.v());
                    flx.z_momentum() += sR*(rho_star*(s_star*nv[2] + qR.w()*(1-nv[2])) - rhoR*qR.w());
                    flx.energy()     += sR*(rho_star*(engyR/rhoR+(s_star-uR_n)*(s_star+qR.p()/(rhoR*(sR-uR_n)))) - engyR);
                }
            }

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
