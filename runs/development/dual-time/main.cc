#include "spade.h"
#include <iomanip>

typedef double real_t;

int main(int argc, char** argv)
{
    real_t q = 1.0;
    real_t rhs = 0.0;
    real_t t0 = 0.0;
    real_t t1 = 5.0*spade::consts::pi;
    int nt = 200;
    real_t dt = (t1 - t0) / (nt);
    q = 5.0;
    real_t t = t0;
    
    struct trans_t
    {
        void transform_forward(real_t& f) const { f = f*f; }
        void transform_inverse(real_t& f) const { f = sqrt(f); }
    } trans;
    
    auto rhs_calc = [&](auto& rhs_in, auto& q_in, const auto& t_in) -> void
    {
        rhs_in = 0.0;
        rhs_in += 2.0*q_in*cos(t_in);
    };
    
    spade::static_math::int_const_t<3> bdf_order;
    
    
    int max_its = 25000;
    const real_t error_tol = 1e-6;
    auto error_norm = [](const real_t& r) -> real_t {return spade::utils::abs(r);};
    spade::time_integration::iterative_control convergence_crit(rhs, error_norm, error_tol, max_its);
    spade::time_integration::dual_time_t time_int(q, rhs, t, dt, dt/10.0, rhs_calc, convergence_crit, bdf_order, trans);
	
    time_int.get_outer_scheme().auxiliary_states[0] = 5.0 + sin(0*dt);
    time_int.get_outer_scheme().auxiliary_states[1] = 5.0 + sin(1*dt);
    time_int.get_outer_scheme().auxiliary_states[2] = 5.0 + sin(2*dt);
    
    trans.transform_forward(time_int.get_outer_scheme().auxiliary_states[0]);
    trans.transform_forward(time_int.get_outer_scheme().auxiliary_states[1]);
    trans.transform_forward(time_int.get_outer_scheme().auxiliary_states[2]);
    
    auto& q0 = time_int.get_outer_scheme().auxiliary_states[0];
    auto& q1 = time_int.get_outer_scheme().auxiliary_states[1];
    auto& q2 = time_int.get_outer_scheme().auxiliary_states[2];
    
    time_int.solution() = 5.0 + sin(2*dt);
    time_int.time() = 2*dt;
    std::ofstream myfile("soln.dat");
    myfile << std::setprecision(15);
    for (auto n: range(0, nt))
    {
        time_int.advance();
        myfile << time_int.time() << " " << time_int.solution() << " " << (5.0+sin(time_int.time())) << std::endl;
    }
    return 0;
}
