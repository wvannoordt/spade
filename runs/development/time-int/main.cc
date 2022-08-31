#include "spade.h"

typedef double real_t;

int main(int argc, char** argv)
{
    real_t rhs, var;
    real_t t0 = 0.0;
    real_t t1 = 5.0*spade::consts::pi;
    int nt = 125;
    real_t dt = (t1 - t0) / (nt);
    var = 5.0;
    auto ftrans = [](real_t& f) -> void {f = f*f; };
    auto itrans = [](real_t& f) -> void {f = sqrt(f); };
    auto calc_rhs = [&](real_t& rhs, const real_t& q, const real_t& time) -> void
    {
        rhs = 0.0;
        rhs += 2.0*q*cos(time);
    };
    spade::time_integration::rk2 time_int(var, rhs, t0, dt, calc_rhs, ftrans, itrans);
    std::ofstream myfile("soln.dat");
    for (int i = 0; i < nt; ++i)
    {
        time_int.advance();
        myfile << time_int.time() << " " << time_int.solution() << " " << (5.0+sin(time_int.time())) << std::endl;
    }
    
    return 0;
}
