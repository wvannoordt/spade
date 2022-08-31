#pragma once

namespace spade::time_integration
{
    template <typename residual_t, typename residual_norm_t, typename error_tol_t> struct iterative_control
    {
        const residual_norm_t* residual_calc;
        error_tol_t error_tol;
        int max_its, current_its;
        typedef error_tol_t err_t;
        
        iterative_control(const residual_t& rhs_in, const residual_norm_t& residual_calc_in, const error_tol_t& error_tol_in, const int& max_its_in)
        {
            max_its = max_its_in;
            error_tol = error_tol_in;
            residual_calc = &residual_calc_in;
            current_its = 0;
        }
        
        error_tol_t calc_residual(const residual_t& r) const { return (*residual_calc)(r); }
    };
}