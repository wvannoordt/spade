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
    
    const real_t fill = 0.0;
    
    // Mesh will tell us how many ODEs to solve
    spade::sym::vector_t vars {u_wm, T_wm, ystar_wm, rho_wm, mu_wm, mut_wm};
    spade::ode::system_t system(fill, spade::ode::constant_mesh_t(y_wm, y, 2), vars, spade::device::best);
    
    using buffer_type = decltype(system)::buffer_type;
    
    spade::ode::expression_t mu_expr(mu_wm, [=] _sp_hybrid (const std::size_t i_inst, const int i_wall, const buffer_type& data)
    {
        const int i = i_wall;
        constexpr real_t T_ref  = 110.4;
        constexpr real_t mu_ref = 1.0e-3;
        const real_t T = data[T_wm][i];
        return mu_ref*pow(T, real_t(1.5))/(T_ref + T); 
    }, spade::ode::xplicit);
    
    const real_t rgas  = 287.15;
    spade::ode::expression_t rho_expr(rho_wm, [=] _sp_hybrid (const std::size_t i_inst, const int i_wall, const buffer_type& data)
    {
        const int i = i_wall;
        const real_t rho_F = data[rho_wm].last();
        return rho_F*(data[T_wm].last()/data[T_wm][i]);
    }, spade::ode::xplicit);
    
    spade::ode::expression_t ystar_expr(ystar_wm, [=] _sp_hybrid (const std::size_t i_inst, const int i_wall, const buffer_type& data)
    {
        // print(i_inst, i_wall);
        const int i = i_wall;
        const real_t mu_w  = data[mu_wm][0];
        const real_t du_dy = (data[u_wm][1] - data[u_wm][0])/(data[y_wm][1] - data[y_wm][0]);
        const real_t tau = mu_w*du_dy;
        return data[y_wm][i]*sqrt(tau*data[rho_wm][0])/data[mu_wm][i];
    }, spade::ode::xplicit);
    
    spade::ode::expression_t mut_expr(mut_wm, [=] _sp_hybrid (const std::size_t i_inst, const int i_wall, const buffer_type& data)
    {
        const int i = i_wall;
        const real_t a_plus = real_t(17.0);
        const real_t kappa  = real_t(0.41);
        const real_t y_coord = data[ystar_wm][i];
        const real_t exp_fac = real_t(1.0) - exp(-y_coord/a_plus);
        return kappa*data[mu_wm][i]*y_coord*exp_fac*exp_fac;
    }, spade::ode::xplicit);
    
    //So on and so forth
    spade::ode::expression_t u_expr (u_wm, [=] _sp_hybrid (const std::size_t i_inst, const int i_wall, const buffer_type& data)
    {
        const int i = i_wall;
        const real_t diff_L = 0.5*(data[mu_wm][i-1] + data[mut_wm][i-1] + data[mu_wm][i] + data[mut_wm][i]);
        const real_t diff_R = 0.5*(data[mu_wm][i+1] + data[mut_wm][i+1] + data[mu_wm][i] + data[mut_wm][i]);
        const real_t dy_L   = data[y_wm][i]   - data[y_wm][i-1];
        const real_t dy_R   = data[y_wm][i+1] - data[y_wm][i];
        const real_t dy_LR  = dy_L + dy_R;
        const real_t c_L    = diff_L/(dy_L*dy_LR);
        const real_t c_R    = diff_R/(dy_R*dy_LR);
        
        real_t lhs0 = c_L;
        real_t lhs1 = -(c_L + c_R);
        real_t lhs2 = c_R;
        real_t rhs  = lhs0*data[u_wm][i-1] + lhs1*data[u_wm][i] + lhs2*data[u_wm][i+1];
        
        return spade::ode::tridiag(lhs0, lhs1, lhs2, rhs);
    }, spade::ode::make_bcs(spade::ode::dirichlet, spade::ode::dirichlet));
    
    const real_t T_F   = 300.0;
    const real_t gamma = 1.4;
    const real_t pr    = 0.71;
    const real_t pr_t  = 0.9;
    spade::ode::expression_t T_expr (T_wm, [=] _sp_hybrid (const std::size_t i_inst, const int i_wall, const buffer_type& data)
    {
        const int i = i_wall;
        const real_t cp     = gamma*rgas/(gamma-real_t(1.0));
        const real_t diff_L = 0.5*(data[mu_wm][i-1]/pr + data[mut_wm][i-1]/pr_t + data[mu_wm][i]/pr + data[mut_wm][i]/pr_t)*cp;
        const real_t diff_R = 0.5*(data[mu_wm][i+1]/pr + data[mut_wm][i+1]/pr_t + data[mu_wm][i]/pr + data[mut_wm][i]/pr_t)*cp;
        const real_t dy_L   = data[y_wm][i]   - data[y_wm][i-1];
        const real_t dy_R   = data[y_wm][i+1] - data[y_wm][i];
        const real_t dy_LR  = dy_L + dy_R;
        const real_t c_L    = diff_L/(dy_L*dy_LR);
        const real_t c_R    = diff_R/(dy_R*dy_LR);
        
        real_t lhs0 = c_L;
        real_t lhs1 = -(c_L + c_R);
        real_t lhs2 = c_R;
        real_t rhs  = lhs0*data[T_wm][i-1] + lhs1*data[T_wm][i] + lhs2*data[T_wm][i+1];
        
        return spade::ode::tridiag(lhs0, lhs1, lhs2, rhs);
    }, spade::ode::make_bcs(spade::ode::zerograd, spade::ode::dirichlet)); //Note that we need to adjust the boundary condition to be e.g. a lambda
    
    auto expressions = spade::ode::make_expressions(rho_expr, mu_expr, ystar_expr, mut_expr, u_expr, T_expr);
    //The order in which these expressions are presented is the order in which they are solved
    
    const real_t uF = 69.54;
    const real_t TF = 300.0;
    system.set_bcs(
        u_wm, [=] _sp_hybrid (const std::size_t& i_inst) { return 69.54; },
        u_wm, [=] _sp_hybrid (const std::size_t& i_inst) { return 0.0;   },
        u_wm, [=] _sp_hybrid (const std::size_t& i_inst) { return 69.54; },
        u_wm, [=] _sp_hybrid (const std::size_t& i_inst) { return 0.0;   },
        );
    
    spade::ode::solve_bvp(system, expressions, spade::utils::converge_crit_t{1.0e-7, 100});
    
    std::ofstream sol("out.dat");
    auto img = system.image();
    for (int i = 0; i < npts; ++i)
    {
        const auto buf = img.get_buffer(0);
        sol << buf[y_wm][i] << " " << buf[u_wm][i] << " " << buf[T_wm][i] << std::endl;
    }
    
    
    return 0;
}
