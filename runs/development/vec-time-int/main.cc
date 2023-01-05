#include "spade.h"

using real_t = double;

int main(int argc, char** argv)
{
    real_t x = 0.1d;
    real_t r = 0.0d;
    
    spade::time_integration::rk4_t alg;
    
    //Note that this data should have the same lifetime as the integrator
    spade::time_integration::integrator_data_t q(x, r, alg);
    
    spade::time_integration::time_axis_t axis(0.0, 1e-4);
    
    //use the "couple" function to build necessary types
    //We can (in theory) extend this to an arbitrary number of integrators
    
    struct trans_t
    {
        void transform_forward(real_t& f) const { f = f*f; }
        void transform_inverse(real_t& f) const { f = sqrt(f); }
    } trans;
    
    auto calc_rhs = [&](auto& rhs, const auto& q, const auto& time) -> void
    {
        rhs = 0.0;
        rhs += 2.0*q*cos(time);
    };
    
    spade::time_integration::integrator_t integrator(
        axis,
        alg,
        q,
        calc_rhs,
        trans
    );
    
    integrator.advance();
    
    // spade::time_integration::rk2 time_int(var, rhs, t0, dt, calc_rhs, trans);
    // std::ofstream myfile("soln.dat");
    // for (int i = 0; i < nt; ++i)
    // {
    //     time_int.advance();
    //     myfile << time_int.time() << " " << time_int.solution() << " " << (5.0+sin(time_int.time())) << std::endl;
    // }
    
    //&spade::ctrs::get<0>(data).solution(0) == &q1.solution(0) //(true!)
    
    return 0;
}
