#pragma once

#include "dispatch/ranges/ranges.h"

namespace spade::dispatch
{
    namespace detail
    {
        template <typename index_t, typename bbx_type, const int idx, typename kernel_t>
        requires ((idx < 0) && ctrs::basic_array<index_t>)
        _sp_hybrid static void cpu_dispatch_impl(index_t& i, const bbx_type& bounds, kernel_t& kernel)
        {
            if constexpr (index_t::size() == 1) kernel(i[0]);
            else                                kernel(i);
        }
        
        template <typename index_t, typename bbx_type, const int idx, typename kernel_t>
        requires ((idx >= 0) && ctrs::basic_array<index_t>)
        _sp_hybrid static void cpu_dispatch_impl(index_t& i, const bbx_type& bounds, kernel_t& kernel)
        {
            for (i[idx] = bounds.min(idx); i[idx] < bounds.max(idx); ++i[idx])
            {
                cpu_dispatch_impl<index_t, bbx_type, idx-1, kernel_t>(i, bounds, kernel);
            }
        }
    }
    
    template <typename int_t, const std::size_t rank, typename idx_t, typename kernel_t>
    _sp_hybrid inline void range_loop(const ranges::basic_range_t<int_t, rank, idx_t>& rg, const kernel_t& kern)
    {
        idx_t i;
        using bound_type = typename utils::remove_all<decltype(rg.bounds)>::type;
        detail::cpu_dispatch_impl<
                idx_t,
                bound_type,
                idx_t::size()-1,
                decltype(kern)>(i, rg.bounds, kern);
    }
    
}