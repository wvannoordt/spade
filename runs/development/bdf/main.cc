#include "spade.h"

typedef double real_t;

int main(int argc, char** argv)
{
    real_t q = 1.0;
    real_t rhs = 0.0;
    real_t t0 = 0.0;
    real_t t1 = 5.0*spade::consts::pi;
    int nt = 500;
    real_t dt = (t1 - t0) / (nt);
    q = 5.0;
    real_t t = t0;
    
    struct trans_t
    {
        void transform_forward(real_t& f) const { f = f*f; }
        void transform_inverse(real_t& f) const { f = sqrt(f); }
    } trans;
    
    
    spade::static_math::int_const_t<3> order;
    auto diff_coeffs = spade::finite_diff::backward_difference_coeffs_node_based<real_t, order.value>();
    const real_t a3 = diff_coeffs[diff_coeffs.size()-1];    
        
    auto rhs_calc = [&](auto& rhs_in, auto& q_in, const auto& t_in) -> void
    {
        rhs_in = 0.0;
        rhs_in += 2.0*q*cos(t_in);
    };
    auto solver   = [&](auto& rhs_in, auto& q_in, auto& rhs_calc_in) -> void
    {
        const real_t k = -dt/a3;
        real_t dr_dq;
        real_t eps = 10.0;
        int it = 0;
        while (spade::utils::abs(eps)>1e-7)
        {
            rhs_calc_in(rhs_in, q_in, t+dt);
            dr_dq = -(2.0*q_in - k*2.0*cos(t+dt));
            q_in = q_in-rhs_in/dr_dq;
            eps = rhs_in;
            it++;
        }
        print(it, eps);
    };
    
    
    spade::time_integration::bdf_t time_int(q, rhs, t, dt, rhs_calc, solver, order, trans);
    time_int.auxiliary_states[0] = 5.0 + sin(0*dt);
    time_int.auxiliary_states[1] = 5.0 + sin(1*dt);
    time_int.auxiliary_states[2] = 5.0 + sin(2*dt);
    trans.transform_forward(time_int.auxiliary_states[0]);
    trans.transform_forward(time_int.auxiliary_states[1]);
    trans.transform_forward(time_int.auxiliary_states[2]);
    q = 5.0 + sin(2*dt);
    trans.transform_forward(q);
    t = 2*dt;
    std::ofstream myfile("soln.dat");
    for (auto n: range(0, nt))
    {
        time_int.advance();
        myfile << time_int.time() << " " << time_int.solution() << " " << (5.0+sin(time_int.time())) << std::endl;
    }
    return 0;
}
