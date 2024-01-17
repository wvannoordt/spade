#pragma once

#include "omni/stencil_union.h"

namespace spade::omni
{
    template <typename... kernels_t>
    struct composite_kernel_t;
    
    template <typename kernel0_t, typename... kernels_t>
    struct composite_kernel_t<kernel0_t, kernels_t...>
    {
        using next_type   = composite_kernel_t<kernels_t...>;
        using output_type = typename kernel0_t::output_type;
        using omni_type   = stencil_union<typename kernel0_t::omni_type, typename next_type::omni_type>;
        kernel0_t kern;
        next_type next;
        
        _sp_hybrid output_type operator() (const auto& input_data) const
        {
            return kern(omni::interpret_stencil<typename kernel0_t::omni_type>(input_data)) + next(omni::interpret_stencil<typename next_type::omni_type>(input_data));
        }
    };
    
    template <typename kernel_t>
    struct composite_kernel_t<kernel_t>
    {
        using output_type = typename kernel_t::output_type;
        using omni_type = typename kernel_t::omni_type;
        
        kernel_t kern;
        
        _sp_hybrid output_type operator() (const auto& input_data) const
        {
            return kern(input_data);
        }
    };
    
    template <typename... kernels_t>
    inline auto compose(const kernels_t&... kernels)
    {
        return composite_kernel_t<kernels_t...>{kernels...};
    }
}