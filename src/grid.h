#pragma once

#include <concepts>
#include <vector>

#include "ctrs.h"
#include "typedef.h"
#include "bounding_box.h"
#include "range.h"
#include "static_for.h"
#include "coord_system.h"
#include "parallel.h"
#include "dims.h"
#include "array_container.h"

namespace cvdf::grid
{
    template <class T> concept multiblock_grid = requires(T t, size_t i, size_t j, size_t k, size_t lb)
    {
        // todo: write this
        { t.node_coords(i, j, k, lb) } -> ctrs::vec_nd<3, typename T::dtype>;
    };
    
    template <coords::coordinate_system coord_t, parallel::parallel_group par_group_t> class cartesian_grid_t
    {
        public:
            typedef coord_t::coord_type dtype;
            cartesian_grid_t(
                const ctrs::array<size_t, cvdf_dim>& num_blocks_in,
                const ctrs::array<size_t, cvdf_dim>& cells_in_block_in,
                const ctrs::array<size_t, cvdf_dim>& exchange_cells_in,
                const bound_box_t<dtype,  cvdf_dim>& bounds_in,
                const coord_t& coord_system_in,
                par_group_t& group_in)
            {
                this->grid_group = &group_in;
                coord_system = coord_system_in;
                dx = 1.0;
                num_blocks = 1;
                cells_in_block = 1;
                total_blocks = 1;
                bounds.min(2) = 0.0;
                bounds.max(2) = 1.0;
                
                for (std::size_t i = 0; i < cvdf_dim; i++)
                {
                    bounds.max(i) = bounds_in.max(i);
                    bounds.min(i) = bounds_in.min(i);
                    num_blocks[i] = num_blocks_in[i];
                    cells_in_block[i] = cells_in_block_in[i];
                    dx[i] = bounds.size(i) / (num_blocks[i]*cells_in_block[i]);
                    exchange_cells[i] = exchange_cells_in[i];
                    total_blocks *= num_blocks[i];
                }
                block_boxes.resize(total_blocks);
                std::size_t clb = 0;
                for (auto lb: range(0,num_blocks[0])*range(0,num_blocks[1])*range(0,num_blocks[2]))
                {
                    auto& box = block_boxes[clb];
                    ctrs::array<dtype, 3> lower;
                    ctrs::array<dtype, 3> upper;
                    static_for<0,3>([&](auto i)
                    {
                        box.min(i.value) = bounds.min(i.value) + (lb[i.value]+0)*bounds.size(i.value)/num_blocks[i.value];
                        box.max(i.value) = bounds.min(i.value) + (lb[i.value]+1)*bounds.size(i.value)/num_blocks[i.value];
                    });
                    ++clb;
                }
            }
            
            ctrs::array<dtype, 3> node_coords(std::size_t i, std::size_t j, std::size_t k, std::size_t lb) const
            {
                ctrs::array<dtype, 3> output(0, 0, 0);
                ctrs::array<std::size_t, 3> ijk(i, j, k);
                auto& box = block_boxes[lb];
                for (size_t idir = 0; idir < cvdf_dim; idir++)
                {
                    output[idir] = box.min(idir) + ijk[idir]*dx[idir];
                }
                return coord_system.map(output);
            }
            
            const par_group_t& group(void) const
            {
                return *grid_group;
            }
            const coord_t& coord_sys(void) const {return coord_system;}
            
            std::size_t get_num_blocks(const std::size_t& i)   const {return num_blocks[i];}
            std::size_t get_num_blocks(void)                   const {return num_blocks[0]*num_blocks[1]*num_blocks[2];}
            std::size_t get_num_cells(const std::size_t& i)    const {return cells_in_block[i];}
            std::size_t get_num_exchange(const std::size_t& i) const {return exchange_cells[i];}
            bound_box_t<dtype,  3> get_bounds(void) const {return bounds;}
            ctrs::array<dtype,  3> get_dx(void) const {return dx;}
            dtype get_dx(const std::size_t& i) const {return dx[i];}
            
        private:
            coord_t coord_system;
            ctrs::array<dtype,  3> dx;
            ctrs::array<std::size_t, 3> num_blocks;
            ctrs::array<std::size_t, 3> cells_in_block;
            ctrs::array<std::size_t, 3> exchange_cells;
            bound_box_t<dtype,  3> bounds;
            std::vector<bound_box_t<dtype, 3>> block_boxes;
            size_t total_blocks;
            par_group_t* grid_group;
    };
    
    template <
        multiblock_grid grid_t,
        typename ar_data_t,
        dims::grid_array_dimension minor_dim_t,
        dims::grid_array_dimension major_dim_t,
        array_container::grid_data_container container_t=std::vector<ar_data_t>>
    struct grid_array
    {
        grid_array(const grid_t& grid_in, const ar_data_t& fill_elem, const minor_dim_t& minor_dims_in, const major_dim_t& major_dims_in)
        {
            minor_dims = minor_dims_in;
            major_dims = major_dims_in;
            grid_dims = dims::dynamic_dims<4>(
                grid_in.get_num_cells(0)+grid_in.get_num_exchange(0)*2,
                grid_in.get_num_cells(1)+grid_in.get_num_exchange(1)*2,
                grid_in.get_num_cells(2)+grid_in.get_num_exchange(2)*2,
                grid_in.get_num_blocks());
            array_container::resize_container(data, minor_dims.total_size()*grid_dims.total_size()*major_dims.total_size());
            array_container::fill_container(data, fill_elem);
            std::size_t n = total_idx_rank();
            for (auto i: range(0,n)) idx_coeffs[i[0]] = 1;
            for (auto i : range(0, n))
            {
                for (auto j : range(i[0]+1, n))
                {
                    idx_coeffs[j[0]] *= this->get_index_extent(i[0]);
                }
            }
            offset = 0;
            for (auto i: range(0,dims::dynamic_dims<4>::rank()-1))
            {
                offset += grid_in.get_num_exchange(i[0])*idx_coeffs[i[0] + minor_dim_t::rank()];
            }
        }
        
        std::size_t get_index_extent(std::size_t i)
        {
            if (i < minor_dim_t::rank())
            {
                return minor_dims.get_index_extent(i);
            }
            else if (i < minor_dim_t::rank()+dims::dynamic_dims<4>::rank())
            {
                return grid_dims.get_index_extent(i-minor_dim_t::rank());
            }
            else
            {
                return major_dims.get_index_extent(i-(minor_dim_t::rank()+dims::dynamic_dims<4>::rank()));
            }
        }
        
        static constexpr std::size_t total_idx_rank(void) {return minor_dim_t::rank()+dims::dynamic_dims<4>::rank()+major_dim_t::rank();}
        
        template <typename idx_t> std::size_t offst_r(const std::size_t& lev, const idx_t& idx) const
        {
            static_assert(std::is_integral<idx_t>::value, "integral type required for array indexing");
            return idx_coeffs[lev]*idx;
        }
        
        template <typename idx_t, typename... idxs_t> std::size_t offst_r(const std::size_t& lev, const idx_t& idx, idxs_t... idxs) const
        {
            static_assert(std::is_integral<idx_t>::value, "integral type required for array indexing");
            return idx_coeffs[lev]*idx + offst_r(lev+1, idxs...);
        }
        
        template <typename... idxs_t>
        ar_data_t& operator() (idxs_t... idxs)
        {
            static_assert(sizeof...(idxs_t)==total_idx_rank(), "wrong number of indices passed to indexing operator");
            return data[offset + offst_r(0, idxs...)];
        }
        
        template <typename... idxs_t>
        const ar_data_t& operator() (idxs_t... idxs) const
        {
            static_assert(sizeof...(idxs_t)==total_idx_rank(), "wrong number of indices passed to indexing operator");
            return data[offset + offst_r(0, idxs...)];
        }
        
        minor_dim_t minor_dims;
        dims::dynamic_dims<4> grid_dims;
        major_dim_t major_dims;
        container_t data;
        std::size_t offset;
        std::size_t idx_coeffs [total_idx_rank()];
    };
}