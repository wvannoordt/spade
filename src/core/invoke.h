#pragma once

#include <concepts>
#include <type_traits>

#include "omni/omni.h"

namespace spade::algs
{
    template <typename kernel_t, typename array_t, typename idx_t>
    auto invoke_at(const array_t& array, const idx_t& i, const kernel_t& kernel)
    {
        //todo: generalize this to a collection of kernels!!

        using omni_type   = kernel_t::omni_type;
        using input_type  = omni::stencil_data_t<omni_type, array_t>;
        input_type data;
        omni::retrieve(array, i, data);
        return kernel(data);
    }
}