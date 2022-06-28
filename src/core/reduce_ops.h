#pragma once

#include <concepts>
#include <iostream>

#include "core/parallel.h"

namespace spade::reduce_ops
{
    template<class T, typename rtype> concept reduce_operation = requires(T t, const rtype b)
    {
        t.value;
        t.init(b);
        t.reduce_elem(b);
        t.equiv_par_op();
    };
    
    template <typename data_t> struct reduce_sum
    {
        data_t value;
        _finline_ void init       (const data_t& init_elem) {value = data_t();}
        _finline_ void reduce_elem(const data_t& new_elem)  {value += new_elem;}
        auto equiv_par_op(void) const {return parallel::par_sum;}
    };
    
    template <typename data_t> struct reduce_max
    {
        data_t value;
        _finline_ void init       (const data_t& init_elem) {value = init_elem;}
        _finline_ void reduce_elem(const data_t& new_elem)  {value = std::max(new_elem,this->value);}
        auto equiv_par_op(void) const {return parallel::par_max;}
    };
}