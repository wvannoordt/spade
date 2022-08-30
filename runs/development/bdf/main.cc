#include "spade.h"

typedef double real_t;

int main(int argc, char** argv)
{
    real_t q = 1.0;
    real_t rhs = 0.0;
    real_t t0 = 0.0;
    real_t dt = 0.01;
    
    auto solver = [](double i) -> double {return 0.0;};
    auto rhs_calc = [](double i) -> double {return 0.0;};
    
    spade::static_math::int_const_t<2> order;
    spade::time_integration::bdf_t time_int(q, rhs, t0, dt, rhs_calc, solver, order);
    for (auto i: range(0,order.value+1))
    {
        print(time_int.diff_coeffs[i]);
    }
    return 0;
}