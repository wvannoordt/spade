#include "spade.h"

int main(int argc, char** argv)
{
    //co-evolving coupled domains
    
    double x1 = 0.1d;
    double r1 = 0.0d;
    
    float  x2 = 0.1f;
    float  r2 = 0.0f;
    
    spade::time_integration::rk2_t alg1;
    spade::time_integration::rk4_t alg2;
    
    auto q1 = spade::time_integration::make_integrator_data(x1, r1, alg1);
    auto q2 = spade::time_integration::make_integrator_data(x2, r2, alg2);
    
    // spade::time_integration::time_axis_t axis(0.0, 1e-4);
    
    // auto rhs1 = [](double& rhs, double& q, const double& t) -> void {rhs = 0.0; rhs -= q;};
    // auto rhs2 = [](float& rhs,  float& q,  const float& t)  -> void {rhs = 0.0; rhs -= q;};
    
    //this is how we would write a single time integrator    
    // spade::time_integration::integrator integrator1(axis, alg1, q1, rhs1);
    // spade::time_integration::integrator integrator2(axis, alg2, q2, rhs2);
    
    //Note that we can't really introduce a composite integrator as a union of two single
    //integrators, since the residual calculation depends on the unions of states
    
    
    // auto coupled_rhs = [](auto& rhs,  auto& q,  const auto& t)  -> void
    // {
    //     auto& resid1 = spade::ctrs::get<0>(rhs);
    //     auto& resid2 = spade::ctrs::get<1>(rhs);
    // 
    //     auto& solut1 = spade::ctrs::get<0>(q);
    //     auto& solut2 = spade::ctrs::get<1>(q);
    // 
    //     resid1 = 0.0;
    //     resid2 = 0.0;
    // 
    //     // etc.
    // };
    
    //use the "couple" function to build necessary types
    //We can (in theory) extend this to an arbitrary number of 
    spade::time_integration::integrator_t c_integrator(
        axis,
        spade::time_integration::couple(alg1, alg2), //const
        spade::time_integration::couple(q1,   q2),   //not const
        coupled_rhs //also need transformation
    );
    
    return 0;
}
