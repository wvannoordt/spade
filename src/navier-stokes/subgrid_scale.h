#pragma once

#include "omni/omni.h"

namespace spade::subgrid_scale
{
    template <typename derived_t, typename ilist_t = omni::info_list_t<>>
    struct subgrid_interface_t
    {
        _sp_hybrid derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        _sp_hybrid const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = ilist_t;

        _sp_hybrid auto get_mu_t(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_mu_t(args...);}, input);
        }
        _sp_hybrid auto get_prt(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_prt(args...);}, input);
        }
    };

    template <typename float_t, typename gas_t> struct wale_t
    : public subgrid_interface_t<wale_t<float_t, gas_t>,
        omni::info_union<omni::info_list_t<omni::info::gradient, omni::info::value>, typename gas_t::info_type>>
    {
        using base_t = subgrid_interface_t<wale_t<float_t, gas_t>,
            omni::info_union<omni::info_list_t<omni::info::gradient, omni::info::value>, typename gas_t::info_type>>;
        using base_t::get_mu_t;
        using base_t::get_prt;
        using base_t::info_type;

        using value_type = float_t;

        float_t cw, delta, prt;
        const gas_t gas;
        wale_t(const gas_t& gas_in, const float_t& cw_in, const float_t& delta_in, const float_t& prt_in)
         : cw{cw_in}, delta{delta_in}, prt{prt_in}, gas{gas_in} {}

        _sp_hybrid float_t get_prt(const auto& info) const {return prt;}

        _sp_hybrid float_t get_mu_t(const auto& info) const
        {
            //for notation, see:
            //Nicoud, F., and Ducros, F.,
            //"Subgrid-Scale Stress Modelling Based on the Square of the Velocity Gradient Tensor,"
            //Flow, Turbulence and Combustion, Vol. 62, No. 3, 1999

            const auto& val       = omni::access<omni::info::value>(info);
            const auto& grad      = omni::access<omni::info::gradient>(info);
            const auto& rgas      = gas.get_R(info);
            const auto  rho       = val.p()/(rgas*val.T());
            const auto  gij       = [&](const int i, const int j){ return grad[j].u(i); };
            const auto  gij2_impl = [&](const int i, const int j){ return gij(i,0)*gij(0,j) + gij(i,1)*gij(1,j) + gij(i,2)*gij(2,j); };
            const auto  sij       = [&](const int i, const int j){ return float_t(0.5)*(gij(i, j) + gij(j, i)); };
            linear_algebra::dense_mat<float_t, 3> gij2(float_t(0.0));
            gij2.fill(gij2_impl);

            linear_algebra::dense_mat<float_t, 3> sdij(float_t(0.0));
            constexpr int dim = sdij.size();
            sdij.fill([&](const int i, const int j){return float_t(0.5)*(gij2(i,j) + gij2(j,i));});

            float_t div = float_t(0.0);
            for (int ii = 0; ii < dim; ++ii) div         += float_t(0.33333333333333333)*gij2(ii,ii);
            for (int ii = 0; ii < dim; ++ii) sdij(ii,ii) -= div;

            const auto tsum = [&](const auto& thing)
            {
                float_t ss = float_t(0.0);
                for (int jj = 0; jj < dim; ++jj)
                {
                    for (int ii = 0; ii < dim; ++ii)
                    {
                        ss += thing(ii,jj)*thing(ii,jj);
                    }
                }
                return ss;
            };

            auto sum_sdij = tsum(sdij);
            auto sum_sij  = tsum(sij);

            const float_t eps   = float_t(1e-8);
            const auto sqrt0    = sqrt(sum_sdij);
            const auto sqrt1    = sqrt(sqrt0);
            const auto sqrt2    = sqrt(sum_sij);
            const auto mu_turb0 = rho*cw*cw*delta*delta*sum_sdij*sqrt0/(eps + sum_sij*sum_sij*sqrt2 + sum_sdij*sqrt1);
            return mu_turb0;
        }
    };
}