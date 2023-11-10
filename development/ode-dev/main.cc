#include "spade.h"

using real_t = double;

int main(int argc, char** argv)
{
    const real_t y0   = 0.0;
    const real_t y1   = 1e-3;
    const real_t dy0  = 1e-6;
    const int    npts = 30;
    const auto   y    = spade::ode::make_stretched_grid(y0, y1, dy0, npts);
    
    using namespace spade::sym::literals;
    
    constexpr static auto u_wm     = "u_wm"_sym;
    constexpr static auto T_wm     = "T_wm"_sym;
    constexpr static auto rho_wm   = "rho_wm"_sym;
    constexpr static auto mu_wm    = "mu_wm"_sym;
    constexpr static auto mut_wm   = "mut_wm"_sym;
    // constexpr static auto zero   = 0_sym;
    
    
    const real_t T_ref  = 110.4;
    const real_t mu_ref = 1.0e-3;
    const auto mu_expression = spade::sym::explcit_t(mu_wm, [=] _sp_hybrid (const int i, const auto& ys, const auto& data)
    {
        const auto& T = data[T_wm][i];
        return mu_ref*pow(T, real_t(1.5))/(T_ref + T);
    });
    
    return 0;
}
