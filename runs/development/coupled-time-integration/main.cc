#include "spade.h"

int main(int argc, char** argv)
{
    double x0 = 1.0d;
    double r0 = 0.0d;
    
    float x1 = 0.5f;
    float r1 = 0.0f;
    
    spade::ctrs::arith_tuple x(x0, x1);
    spade::ctrs::arith_tuple r(r0, r1);

    // spade::time_integration::ssprk3_t alg;
    // spade::time_integration::rk38r_t alg;
    // spade::time_integration::ssprk34_t alg;
    // spade::time_integration::rk4_t alg;
    spade::time_integration::rk2_t alg;
    
    //Note that this data should have the same lifetime as the integrator
    spade::time_integration::integrator_data_t q(x, r, alg);
    
    spade::time_integration::time_axis_t axis(0.0, 1e-2);

    auto calc_rhs = [&](auto& rhs, const auto& sol, const auto& time) -> void
    {
        auto& rhs0 = spade::ctrs::get<0>(rhs);
        auto& rhs1 = spade::ctrs::get<1>(rhs);
        
        auto& q0 = spade::ctrs::get<0>(sol);
        auto& q1 = spade::ctrs::get<1>(sol);
        
        rhs0 = 0.0;
        rhs1 = 0.0;

        rhs0 += 0.01*q0 - 0.05*q1;
        rhs1 += 0.04*q0 + 0.03*q1;
    };
    
    spade::time_integration::integrator_t integrator(
        axis,
        alg,
        q,
        calc_rhs
    );
    
    integrator.advance();
    
    std::ofstream myfile("soln.dat");
    myfile << std::setprecision(15);
    const int nt = 20000;
    for (int i = 0; i < nt; ++i)
    {
        const auto& soln = integrator.solution();
        const auto& numsol0  = spade::ctrs::get<0>(soln);
        const auto& numsol1  = spade::ctrs::get<1>(soln);
        myfile
            << integrator.time() << " "
            << numsol0 << " "
            << numsol1 << "\n";
        integrator.advance();
    }
    
    return 0;
}
