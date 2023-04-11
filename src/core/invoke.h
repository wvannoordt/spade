#pragma once

#include <concepts>
#include <type_traits>

#include "omni/omni.h"

namespace spade::algs
{
    template <typename array_t, typename idx_t, typename kernel_t>
    auto invoke_at(const array_t& array, const idx_t& i, const kernel_t& kernel)
    {
        //todo: generalize this to a collection of kernels!!
        using in_type          = kernel_t::omni_type;
        using adj_type         = in_type::template recenter<idx_t::centering_type()>;
        constexpr bool req_ctr = in_type::center() == grid::agno_centered;
        using omni_type   = typename std::conditional<req_ctr, adj_type, in_type>::type;
        using input_type  = omni::stencil_data_t<omni_type, array_t>;
        input_type data;
        omni::retrieve(array, i, data);
        return kernel(data);
    }
}