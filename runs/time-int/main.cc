#include "cvdf.h"

typedef double real_t;

int main(int argc, char** argv)
{
    real_t rhs1, rhs2, var;
    real_t t0 = 0.0;
    real_t t1 = 5.0*cvdf::consts::pi;
    int nt = 250;
    real_t dt = (t1 - t0) / (nt);
    var = -1.0;
    
    auto calc_rhs = [&](real_t& rhs, const real_t& q, const real_t& time) -> void
    {
        rhs = 0.0;
        rhs += sin(time);
    };
    cvdf::time_integration::identity_transform_t<real_t> trans;
    cvdf::time_integration::rk2 time_int(var, rhs1, rhs2, t0, dt, calc_rhs, trans, trans);
    std::ofstream myfile("soln.dat");
    for (int i = 0; i < nt; ++i)
    {
        time_int.advance();
        print(time_int.time(), time_int.soln(), -cos(time_int.time()));
        myfile << time_int.time() << " " << time_int.soln() << " " << -cos(time_int.time()) << std::endl;
    }
    
    return 0;
}