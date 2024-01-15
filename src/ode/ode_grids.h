#pragma once

#include "core/num_algs.h"
#include "dispatch/shared_vector.h"

namespace spade::ode
{
    template <typename float_t>
    device::shared_vector<float_t> make_stretched_grid(const float_t& x0, const float_t& x1, const float_t& dx, const int npts)
    {
        device::shared_vector<float_t> output;
        output.resize(npts);
        const int n = npts - 1;
        const float_t xn = x1 - x0;
        
        
        constexpr int max_its = 100;
        constexpr float_t tol = float_t(1e-9);
        const auto fx   = [&](const float_t& r) { return float_t(1.0) - std::pow(r, n) - (xn/dx)*(float_t(1.0) - r); };
        const auto dfdx = [&](const float_t& r) { return - n*std::pow(r, n-1) + (xn/dx); };
        const auto result = num_algs::newton(float_t(2.0), fx, dfdx, max_its, tol);
        auto ratio = result.x;
        auto its   = result.its;
        auto err   = result.eps;
        if (its > max_its)
        {
            std::stringstream msg;
            msg << "Failed ODE stretched grid generation for:\nx0 = " << x0 << "\n";
            msg << "x1 = " << x1 << "\ndx = " << dx << "\nn  = " << npts;
            throw except::sp_exception(msg.str());
        }
        output[0] = x0;
        float_t dx_loc = dx;
        for (int i = 1; i < n; ++i)
        {
            output[i] = output[i-1] + dx_loc;
            dx_loc *= ratio;
        }
        output[n] = x1;
        output.transfer();
        return output;
    }
}