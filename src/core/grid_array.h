#pragma once

#include <concepts>
#include <vector>
#include <type_traits>

#include "core/attribs.h"
#include "core/ctrs.h"
#include "core/typedef.h"
#include "core/bounding_box.h"
#include "core/range.h"
#include "core/static_for.h"
#include "core/coord_system.h"
#include "core/parallel.h"
#include "core/dims.h"
#include "core/array_container.h"
#include "core/partition.h"
#include "core/static_math.h"
#include "core/grid_index_types.h"
#include "core/grid.h"

namespace spade::grid
{    
    namespace detail
    {
        template <class T> concept has_value_type = requires(T t) {typename T::value_type;};
        
        template <typename data_t> struct get_fundamental_type
        {
            typedef data_t type;
        };
        
        template <has_value_type data_t> struct get_fundamental_type<data_t>
        {
            typedef typename data_t::value_type type;
        };
        
        template <class T> concept is_static_1D_array = ctrs::basic_array<T> && requires (T t)
        {
            t; //todo: figure out whay the heck is going on with this
        };        
        template <typename data_t> struct get_dim_type
        {
            typedef dims::singleton_dim type;
        };
        
        template <is_static_1D_array data_t> struct get_dim_type<data_t>
        {
            typedef dims::static_dims<data_t::size()> type;
        };
    }
    
    template <
        multiblock_grid grid_t,
        typename data_alias_t,
        dims::grid_array_dimension major_dim_t=dims::singleton_dim,
        const array_center_e centering=cell_centered,
        array_container::grid_data_container container_t=std::vector<typename detail::get_fundamental_type<data_alias_t>::type>>
    struct grid_array
    {
        typedef typename detail::get_fundamental_type<data_alias_t>::type fundamental_type;
        typedef typename detail::get_dim_type<data_alias_t>::type minor_dim_t;
        grid_array(){}
        grid_array(
            const grid_t& grid_in,
            const data_alias_t& fill_elem,
            const major_dim_t& major_dims_in=dims::singleton_dim())
        {
            grid = &grid_in;
            minor_dims = minor_dim_t();
            major_dims = major_dims_in;
            grid_dims = create_grid_dims(grid_in, this->centering_type());
            array_container::resize_container(data, minor_dims.total_size()*grid_dims.total_size()*major_dims.total_size());
            array_container::fill_container(data, fill_elem);
            std::size_t n = total_idx_rank();
            for (auto i: range(0,n)) idx_coeffs[i] = 1;
            for (auto i : range(0, n))
            {
                for (int j =i+1; j < n; ++j)
                {
                    idx_coeffs[j] *= this->get_index_extent(i);
                }
            }
            offset = 0;
            for (auto i: range(0,dims::dynamic_dims<4>::rank()-1))
            {
                offset += grid_in.get_num_exchange(i)*idx_coeffs[i + minor_dim_t::rank()];
            }
        }
        
        typedef grid_t grid_type;
        typedef fundamental_type value_type;
        typedef minor_dim_t array_minor_dim_t;
        typedef major_dim_t array_major_dim_t;
        typedef std::conditional<centering==cell_centered, cell_t<int>, node_t<int>> index_integral_t;
        typedef data_alias_t alias_type;
        typedef data_alias_t unwrapped_minor_type;
        
        constexpr static array_center_e centering_type(void) { return centering; }
        
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
        
        template <typename idx_t>
        requires std::convertible_to<idx_t, int>
         _finline_ std::size_t offst_r(const std::size_t& lev, const idx_t& idx) const
        {
            return idx_coeffs[lev]*idx;
        }
        
        template <typename idx_t, typename... idxs_t>
        requires std::convertible_to<idx_t, int>
        _finline_ std::size_t offst_r(const std::size_t& lev, const idx_t& idx, idxs_t... idxs) const
        {
            return idx_coeffs[lev]*idx + offst_r(lev+1, idxs...);
        }
        
        template <typename... idxs_t>
        _finline_ fundamental_type& operator() (idxs_t... idxs)
        {
            static_assert(sizeof...(idxs_t)==total_idx_rank(), "wrong number of indices passed to indexing operator");
            return data[offset + offst_r(0, idxs...)];
        }
        
        template <typename... idxs_t>
        _finline_ const fundamental_type& operator() (idxs_t... idxs) const
        {
            static_assert(sizeof...(idxs_t)==total_idx_rank(), "wrong number of indices passed to indexing operator");
            return data[offset + offst_r(0, idxs...)];
        }
        
        template <typename i0_t, typename i1_t, typename i2_t, typename i3_t, typename i4_t, typename i5_t>
        requires std::convertible_to<i0_t, int> &&
        std::convertible_to<i1_t, int> &&
        std::convertible_to<i2_t, int> &&
        std::convertible_to<i3_t, int> &&
        std::convertible_to<i4_t, int> &&
        std::convertible_to<i5_t, int>
        _finline_ fundamental_type& unwrap_idx(const i0_t& i0, const i1_t& i1, const i2_t& i2, const i3_t& i3, const i4_t& i4, const i5_t& i5)
        {
            return data[
                offset +
                i0*idx_coeffs[0] + 
                i1*idx_coeffs[minor_dim_t::rank()+0]+
                i2*idx_coeffs[minor_dim_t::rank()+1]+
                i3*idx_coeffs[minor_dim_t::rank()+2]+
                i4*idx_coeffs[minor_dim_t::rank()+3]+
                i5*idx_coeffs[minor_dim_t::rank()+4]
            ];
        }
        
        template <typename i0_t, typename i1_t, typename i2_t, typename i3_t, typename i4_t, typename i5_t>
        _finline_ const fundamental_type& unwrap_idx(const i0_t& i0, const i1_t& i1, const i2_t& i2, const i3_t& i3, const i4_t& i4, const i5_t& i5) const
        {
            static_assert(std::is_integral<i0_t>::value, "unwrap_idx requires integral arguments");
            static_assert(std::is_integral<i1_t>::value, "unwrap_idx requires integral arguments");
            static_assert(std::is_integral<i2_t>::value, "unwrap_idx requires integral arguments");
            static_assert(std::is_integral<i3_t>::value, "unwrap_idx requires integral arguments");
            static_assert(std::is_integral<i4_t>::value, "unwrap_idx requires integral arguments");
            static_assert(std::is_integral<i5_t>::value, "unwrap_idx requires integral arguments");
            return data[
                offset +
                i0*idx_coeffs[0] + 
                i1*idx_coeffs[minor_dim_t::rank()+0]+
                i2*idx_coeffs[minor_dim_t::rank()+1]+
                i3*idx_coeffs[minor_dim_t::rank()+2]+
                i4*idx_coeffs[minor_dim_t::rank()+3]+
                i5*idx_coeffs[minor_dim_t::rank()+4]
            ];
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator -= (const rhs_t& rhs)
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] -= rhs.data[i];
            return *this;
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator += (const rhs_t& rhs)
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] += rhs.data[i];
            return *this;
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator *= (const rhs_t& rhs)
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] *= rhs.data[i];
            return *this;
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator /= (const rhs_t& rhs)
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] /= rhs.data[i];
            return *this;
        }
        
        template <typename numeric_t> grid_array& operator *= (const numeric_t& rhs)
        requires std::floating_point<numeric_t>
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] *= rhs;
            return *this;
        }
        
        template <typename numeric_t> grid_array& operator /= (const numeric_t& rhs)
        requires std::floating_point<numeric_t>
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] /= rhs;
            return *this;
        }
        
        grid_array& operator /= (const grid_array& rhs)
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] /= rhs.data[i];
            return *this;
        }
        
        grid_array& operator = (const fundamental_type& v)
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] = v;
            return *this;
        }
        
        minor_dim_t get_minor_dims(void) const { return minor_dims; }
        major_dim_t get_major_dims(void) const { return major_dims; }
        major_dim_t get_grid_dims (void) const { return grid_dims;  }
        auto get_total_dims(void)        const { return minor_dims*grid_dims*major_dims; }
        const grid_t& get_grid(void) const {return *grid;}
        
        const grid_t* grid;
        minor_dim_t minor_dims;
        dims::dynamic_dims<4> grid_dims;
        major_dim_t major_dims;
        container_t data;
        std::size_t offset;
        std::size_t idx_coeffs [total_idx_rank()];
    };
}