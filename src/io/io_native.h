#pragma once

#include <iostream>
#include <fstream>

#include "core/config.h"
#include "core/static_for.h"
#include "grid/grid.h"
#include "core/print.h"
#include "core/utils.h"
#include "core/base_64.h"
#include "core/mem_map.h"

namespace spade::io
{
    namespace detail
    {
        template <grid::multiblock_array array_t> void get_array_par_buf(parallel::par_buf_t& buf, const array_t& array)
        {
            buf.clear();
            const auto& grid = array.get_grid();
            std::size_t block_elems = mem_map::map_size(array.var_map())*mem_map::map_size(array.block_map()) / grid.get_num_local_blocks();
            using data_t = typename array_t::value_type;
            std::size_t block_size_bytes = block_elems*sizeof(data_t);
            ctrs::array<int, 3> nexch (grid.get_num_exchange(0), grid.get_num_exchange(1), grid.get_num_exchange(2));
            for (auto lb: range(0, grid.get_num_local_blocks()))
            {
                std::size_t lb_glob = array.get_grid().get_partition().get_global_block(lb);
                ctrs::array<int, array_t::variable_map_type::rank()> iv = 0;
                void* ptr = (void*)(&array(iv,-nexch[0], -nexch[1], -nexch[2], lb));
                std::size_t offset_bytes = lb_glob*block_size_bytes;
                buf.add(ptr, block_size_bytes, offset_bytes);
            }
        }
    }
    
    template <grid::multiblock_array array_t> void binary_write(const std::string& filename, array_t& array)
    {
        const auto& group = array.get_grid().group();
        parallel::par_buf_t buf;
        detail::get_array_par_buf(buf, array);
        parallel::mpi_file_t mf(group, filename);
        mf.write_buf(buf);
    }
    
    template <grid::multiblock_array array_t> void binary_read(const std::string& filename, array_t& array)
    {
        const auto& group = array.get_grid().group();
        parallel::par_buf_t buf;
        detail::get_array_par_buf(buf, array);
        parallel::mpi_file_t mf(group, filename);
        mf.read_buf(buf);
    }
}