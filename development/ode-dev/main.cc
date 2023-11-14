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
    
    constexpr static auto y_wm     = "y_wm"_sym;
    constexpr static auto ystar_wm = "ystar_wm"_sym;
    constexpr static auto u_wm     = "u_wm"_sym;
    constexpr static auto T_wm     = "T_wm"_sym;
    constexpr static auto rho_wm   = "rho_wm"_sym;
    constexpr static auto mu_wm    = "mu_wm"_sym;
    constexpr static auto mut_wm   = "mut_wm"_sym;
    // constexpr static auto zero   = 0_sym; // This is nice!
    
    // CTAD for underlying value type (fundamental only)
    const real_t fill = 0.0;
    
    // Mesh will tell us how many ODEs to solve
    spade::sym::vector_t vars {u_wm, T_wm, ystar_wm, rho_wm, mu_wm, mut_wm};
    // spade::ode::system_t system(fill, spade::ode::constant_mesh_t(y_wm, y, 100), vars);
    /*
    using buffer_type = decltype(system)::buffer_type;
    
    spade::ode::expression_t mu_expr(mu_wm, [=] _sp_hybrid (const int i, const buffer_type& data)
    {
        constexpr real_t T_ref  = 110.4;
        constexpr real_t mu_ref = 1.0e-3;
        const auto& T = data[T_wm][i];
        return mu_ref*pow(T, real_t(1.5))/(T_ref + T);
    });
    
    spade::ode::expression_t rho_expr(rho_wm, [=] _sp_hybrid (const int i, const buffer_type& data)
    {
        return data[rho_wm][npts];
    });
    
    spade::ode::expression_t ystar_expr(ystar_wm, [=] _sp_hybrid (const int i, const buffer_type& data)
    {
        const real_t mu_w  = data[mu_wm][0];
        const real_t du_dy = (data[u_wm][1] - data[u_wm][0])/(data[y_wm][1] - data[y_wm][0]);
        const real_t tau = mu_w*du_dy;
        return data[y_wm][]*sqrt(tau*data[rho_wm][0])/data[mu_wm][i];
    });
    
    spade::ode::expression_t mut_expr(mut_wm, [=] _sp_hybrid (const int i, const buffer_type& data)
    {
        const real_t a_plus = real_t(17.0);
        const real_t kappa  = real_t(0.41);
        const real_t y_coord = data[ystar_wm][i];
        const real_t exp_fac = real_t(1.0) - exp(-y_coord/a_plus);
        return kappa*data[mu_wm][i]*y_coord*exp_fac;
    });
    
    //So on and so forth
    spade::ode::expression_t u_expr (u_wm, [=] _sp_hybrid (const int i, const buffer_type& data)
    {
        return 0.0;
    });
    
    spade::ode::expression_t T_expr (T_wm, [=] _sp_hybrid (const int i, const buffer_type& data)
    {
        return 0.0;
    });
    
    //Boundary conditions (tbd)
    // system.set_boundary([=] _sp_hybrid(const std::size_t i){ return 0.0; });
    
    //Initialize data
    // system.initialize(u_wm, []);
    
    //The order in which these expressions are presented is the order in which they are solved
    
    spade::ode::solve_bvp(system, {rho_expr, ystar_expr, mu_expr, mut_expr, T_expr, u_expr});
    */
    return 0;
}
