#pragma once

namespace spade::time_integration
{
    //represents the fixed point iteration:
    //
    //  w_{j+1} = coeff*w_{j} + (1-coeff)*(w_0 + 0.5*dt*R(w_0) + R(w_{j}))
    //
    template <typename err_func_t, typename coeff_t>
    struct crank_nicholson_t
    {
        using coeff_type = coeff_t;
        coeff_t coeff;
        err_func_t err_func;
        utils::converge_crit_t<coeff_t> crit;
        crank_nicholson_t(const err_func_t& err_func_in, const utils::converge_crit_t<coeff_t>& crit_in, const coeff_t& coeff_in = 0.875)
        : err_func{err_func_in}, crit{crit_in}, coeff{coeff_in} {}

        static constexpr std::size_t var_size() {return 1;}
        static constexpr std::size_t rhs_size() {return 2;}

        constexpr static bool is_crank_nichol_specialization = true;
    };
}
