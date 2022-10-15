#pragma once

namespace spade::grid
{
    template <typename something_t> struct md_mem_map
    {
        typedef int index_type;
        typedef std::size_t offset_type;
        
        template <typename... idxs_t> offset_type compute_offset(const idxs_t&... idxs)
        {
            return 0;
        }
    }
}