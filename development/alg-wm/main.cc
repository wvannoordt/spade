#include "spade.h"
#include "dispatch/execute2.h"
using real_t = double;

int main(int argc, char** argv)
{
    const auto device = spade::device::best;
    const int n = 5;
    spade::device::shared_vector<real_t> u_in;
    spade::device::shared_vector<real_t> rho_in;
    spade::device::shared_vector<real_t> y_in;
    u_in.resize(n);
    rho_in.resize(n);
    y_in.resize(n);
    
    for (auto& u: u_in.data  (spade::device::cpu)) u = 69.54 + 3.5*spade::utils::unitary_random();
    for (auto& r: rho_in.data(spade::device::cpu)) r = 1.17  + 0.0056*spade::utils::unitary_random();
    for (auto& y: y_in.data  (spade::device::cpu)) y = 0.001  + 0.0001*spade::utils::unitary_random();
    
    u_in.transfer();
    rho_in.transfer();
    y_in.transfer();
    
    const real_t mu    = 1.8e-5;
    const auto u_img   = spade::utils::make_vec_image(u_in.data(device));
    const auto rho_img = spade::utils::make_vec_image(rho_in.data(device));
    const auto y_img   = spade::utils::make_vec_image(y_in.data(device));
    
    spade::device::shared_vector<real_t> tau_out;
    tau_out.resize(n);
    
    auto tau_img = spade::utils::make_vec_image(tau_out.data(device));
    
    auto loop = [=] _sp_hybrid (const std::size_t i) mutable
    {
        const real_t rho_loc = rho_img[i];
        const real_t y_loc   = y_img[i];
        const real_t u_loc   = u_img[i];
        
        const auto func = [&](const real_t& u_tau)
        {
            const auto upls = [&](const real_t& y_plus)
            {
                constexpr real_t c_b  = 5.0333908790505579;
                constexpr real_t c_a1 = 8.148221580024245;
                constexpr real_t c_a2 = -6.9287093849022945;
                constexpr real_t c_b1 = 7.4600876082527945;
                constexpr real_t c_b2 = 7.468145790401841;
                constexpr real_t c_c1 = 2.5496773539754747;
                constexpr real_t c_c2 = 1.3301651588535228;
                constexpr real_t c_c3 = 3.599459109332379;
                constexpr real_t c_c4 = 3.6397531868684494;
                
                return c_b + c_c1*log((y_plus+c_a1)*(y_plus+c_a1)+c_b1*c_b1) - c_c2*log((y_plus+c_a2)*(y_plus+c_a2)+c_b2*c_b2)
	               -c_c3*atan2(c_b1, y_plus+c_a1)-c_c4*atan2(c_b2, y_plus+c_a2);
            };
            
            const real_t u_plus_end0 = upls(rho_loc*u_tau*y_loc/mu);
            const real_t u_plus_end1 = u_loc/u_tau;
            
            return u_plus_end0 - u_plus_end1;
        };
        
        
        
        const real_t u_tau_ini = u_img[i]*0.1;
        const auto result = spade::num_algs::newton(u_tau_ini, func, spade::num_algs::num_deriv(1.0e-6), 100, real_t(1.0e-6));
        tau_img[i] = rho_loc*result.x*result.x;
    };
    
    //This should no longer take a device
    auto range = spade::dispatch::ranges::from_array(tau_img, device);
    spade::dispatch::execute(range, loop);
    
    tau_out.itransfer();
    
    for (auto tau: tau_out)
    {
        print(tau);
    }
    
    return 0;
}
