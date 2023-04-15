#pragma once

#include "omni/omni.h"

namespace spade::subgrid_scale
{
    template <typename derived_t, omni::info::is_omni_info... infos_t>
    struct subgrid_interface_t
    {
        derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = omni::info_list_t<infos_t...>;

        auto get_mut(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_nu_turb(args...);}, input);
        }
    };

    template <typename float_t> struct wale_t
    : public subgrid_interface_t<wale_t<float_t>, omni::info::gradient>
    {
        using base_t = subgrid_interface_t<wale_t<float_t>, omni::info::gradient>;
        using base_t::get_mut;
        using base_t::info_type;

        float_t cw, delta;
        wale_t(const float_t& cw_in, const float_t& delta_in) : cw{cw_in}, delta{delta_in} {}

        float_t get_nu_turb(const ctrs::array<fluid_state::prim_t<float_t>, 3>& grad) const
        {
            //for notation, see:
            //Nicoud, F., and Ducros, F.,
            //"Subgrid-Scale Stress Modelling Based on the Square of the Velocity Gradient Tensor,"
            //Flow, Turbulence and Combustion, Vol. 62, No. 3, 1999

            const auto gij       = [&](const int i, const int j){ return grad[j].u(i); };
            const auto gij2_impl = [&](const int i, const int j){ return gij(i,0)*gij(0,j) + gij(i,1)*gij(1,j) + gij(i,2)*gij(2,j); };
            const auto sij       = [&](const int i, const int j){ return float_t(0.5)*(gij(i, j) + gij(j, i)); };
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

            const float_t eps = float_t(1e-8);
            auto nu_turb = cw*cw*delta*delta*std::pow(sum_sdij, 1.5)/(eps + std::pow(sum_sij, 2.5) + std::pow(sum_sdij, 1.25));
            return nu_turb;
        }
    };
}