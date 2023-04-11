#pragma once

#include "core/typedef.h"

namespace spade::utils
{
    template <typename data_t, typename coeff_t = real_t>
    struct avg_t
    {
        data_t val, val2;
        std::size_t count;
        avg_t() : val{data_t()}, val2{data_t()*data_t()} {}
        avg_t& operator << (const data_t& nval)
        {
            const coeff_t alph = coeff_t(1.0)/(count+1);
            const coeff_t beta = coeff_t(1.0)-alph;
            val  =  alph*nval      + beta*val;
            val2 = alph*nval*nval + beta*val2;
            ++count;
            return *this;
        }

        data_t value() const {return val;}
        data_t var()   const {return val2-val*val;}
        data_t stdev() const {return std::sqrt(var());}
    };
}