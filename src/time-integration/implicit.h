#pragma once

namespace spade::time_integration
{
    //represents the fixed point iteration:
    //
    //  w_{j+1} = coeff*w_{j} + (1-coeff)*(w_0 + 0.5*dt*R(w_0) + R(w_{j}))
    //
    template <typename coeff_t>
    struct crank_nicholson_t
    {
        using coeff_type = coeff_t;
        coeff_t coeff;
        crank_nicholson_t() : coeff{0.875} {}
        crank_nicholson_t(const coeff_t& coeff_in) : coeff{coeff_in} {}

        static constexpr std::size_t var_size() {return 2;}
        static constexpr std::size_t rhs_size() {return 2;}

        constexpr static bool is_crank_nichol_specialization = true;
    };
}