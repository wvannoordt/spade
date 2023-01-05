#include "spade.h"

using real_t = double;

int main(int argc, char** argv)
{
    real_t x = 5.0d;
    real_t r = 0.0d;
    
    spade::time_integration::ssprk34_t alg;
    
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
    
    std::ofstream myfile("soln.dat");
    myfile << std::setprecision(15);
    const int nt = 100000;
    for (int i = 0; i < nt; ++i)
    {
        integrator.advance();
        myfile << integrator.time() << " " << integrator.solution() << " " << (5.0+sin(integrator.time())) << std::endl;
    }
    
    return 0;
}
