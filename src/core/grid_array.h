#pragma once

#include <vector>
#include <type_traits>

#include "core/config.h"
#include "core/attribs.h"
#include "core/ctrs.h"
#include "core/typedef.h"
#include "core/bounding_box.h"
#include "core/range.h"
#include "core/coord_system.h"
#include "core/parallel.h"
#include "core/dims.h"
#include "core/partition.h"
#include "core/static_math.h"
#include "core/grid_index_types.h"
#include "core/grid.h"
#include "core/mem_map.h"

// #include "core/array_container.h" //deprecated
#include "array-containers/ac_vector_wrapper.h"

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
                
        //todo: phase this out
        template <typename data_t> struct get_dim_type
        {
            typedef dims::singleton_dim type;
        };
        
        template <is_static_1D_array data_t> struct get_dim_type<data_t>
        {
            typedef dims::static_dims<data_t::size()> type;
        };
        
        template <typename data_t> struct get_variable_mem_map
        {
            typedef singleton_map_t type;
        };
        
        template <is_static_1D_array data_t> struct get_variable_mem_map<data_t>
        {
            typedef regular_map_t<static_dim_t<0,data_t::size()>> type;
        };
        
        
        //getting the memory map type for the grid indices
        template <const array_centering center, const std::size_t grid_dim> struct get_ijklb_map_type{};
        
        template <const std::size_t grid_dim> struct get_ijklb_map_type<cell_centered, grid_dim>
        {
            typedef regular_map_t<
                dynamic_dim_t<int>,
                dynamic_dim_t<int>,
                dynamic_dim_t<int>,
                dynamic_dim_t<int>> type;
        };
        
        template <const std::size_t grid_dim> struct get_ijklb_map_type<node_centered, grid_dim>
        {
            typedef regular_map_t<
                dynamic_dim_t<int>,
                dynamic_dim_t<int>,
                dynamic_dim_t<int>,
                dynamic_dim_t<int>> type;
        };
        
        template <const std::size_t grid_dim> struct get_ijklb_map_type<face_centered, grid_dim>
        {
            //Can I do a nicer job here?
            typedef regular_map_t<
                typename std::conditional<0==face_idx_t::dir_idx, static_dim_t<0, grid_dim>, dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, static_dim_t<0, grid_dim>, dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, static_dim_t<0, grid_dim>, dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, static_dim_t<0, grid_dim>, dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, static_dim_t<0, grid_dim>, dynamic_dim_t<int>>::type
                > type;
        };
        
        template <const array_centering centering> struct get_centering_grid_idx_rank
        {
            static constexpr int value = 4;
        };
        
        template <> struct get_centering_grid_idx_rank<cell_centered> { static constexpr int value = 4; };
        template <> struct get_centering_grid_idx_rank<face_centered> { static constexpr int value = 5; };
        template <> struct get_centering_grid_idx_rank<node_centered> { static constexpr int value = 4; };
        
        template <const array_centering centering> struct get_index_type{};
        template <> struct get_index_type<cell_centered>
        {
            typedef typename cell_idx_t::value_type integral_type;
            typedef cell_idx_t array_type;
        };
        template <> struct get_index_type<face_centered>
        {
            typedef typename face_idx_t::value_type integral_type;
            typedef face_idx_t array_type;
        };
        template <> struct get_index_type<node_centered>
        {
            typedef typename node_idx_t::value_type integral_type;
            typedef node_idx_t array_type;
        };
        
        template <typename idx_t> 
        struct count_rank_t
        {
            static constexpr std::size_t value = 1;
        };
        
        template <ctrs::basic_array idx_t> 
        struct count_rank_t<idx_t>
        {
            static constexpr std::size_t value = idx_t::size();
        };
        
        template <typename idx_t, typename... idxs_t> struct index_rank_size_t
        {
            static constexpr std::size_t value = count_rank_t<idx_t>::value + index_rank_size_t<idxs_t...>::value;
        };
        
        template <typename idx_t> struct index_rank_size_t<idx_t>
        {
            static constexpr std::size_t value = count_rank_t<idx_t>::value;
        };
        
        template <const array_centering centering, typename grid_map_t, typename grid_t>
        requires (centering == cell_centered)
        void insert_grid_dims(grid_map_t& map, const grid_t& grid)
        {
            //i
            std::get<cell_idx_t::i_idx>(map.dims).i0 = -grid.get_num_exchange(0);
            std::get<cell_idx_t::i_idx>(map.dims).i1 =  grid.get_num_cells(0) + grid.get_num_exchange(0);
            
            //j
            std::get<cell_idx_t::j_idx>(map.dims).i0 = -grid.get_num_exchange(1);
            std::get<cell_idx_t::j_idx>(map.dims).i1 =  grid.get_num_cells(1) + grid.get_num_exchange(1);
            
            //k
            std::get<cell_idx_t::k_idx>(map.dims).i0 = -grid.get_num_exchange(2);
            std::get<cell_idx_t::k_idx>(map.dims).i1 =  grid.get_num_cells(2) + grid.get_num_exchange(2);
            
            //lb
            std::get<cell_idx_t::lb_idx>(map.dims).i0 = 0;
            std::get<cell_idx_t::lb_idx>(map.dims).i1 = grid.get_num_local_blocks();
            
            map.compute_coeffs();
        }
        
        template <const array_centering centering, typename grid_map_t, typename grid_t>
        requires (centering == face_centered)
        void insert_grid_dims(grid_map_t& map, const grid_t& grid)
        {
            //i
            std::get<face_idx_t::i_idx>(map.dims).i0 = -grid.get_num_exchange(0);
            std::get<face_idx_t::i_idx>(map.dims).i1 =  grid.get_num_cells(0) + grid.get_num_exchange(0) + 1;
            
            //j
            std::get<face_idx_t::j_idx>(map.dims).i0 = -grid.get_num_exchange(1);
            std::get<face_idx_t::j_idx>(map.dims).i1 =  grid.get_num_cells(1) + grid.get_num_exchange(1) + 1;
            
            //k
            std::get<face_idx_t::k_idx>(map.dims).i0 = -grid.get_num_exchange(2);
            std::get<face_idx_t::k_idx>(map.dims).i1 =  grid.get_num_cells(2) + grid.get_num_exchange(2) + 1;
            
            //lb
            std::get<face_idx_t::lb_idx>(map.dims).i0 = 0;
            std::get<face_idx_t::lb_idx>(map.dims).i1 = grid.get_num_local_blocks();
            
            map.compute_coeffs();
        }
        
    }
    
    namespace detail
    {
        template <const array_centering centering, typename bbox_t>
        requires (centering == face_centered)
        void set_bbox_dir_val(bbox_t& bbox, const int& dim)
        {
            bbox.min(face_idx_t::dir_idx) = 0;
            bbox.max(face_idx_t::dir_idx) = dim;
        }
        
        template <const array_centering centering, typename bbox_t>
        requires (centering != face_centered)
        void set_bbox_dir_val(bbox_t& bbox, const int& dim) {}
    }
    
    template <typename data_t> using default_container_t = array_containers::vector_wrapper<data_t>;
    
    template <
        multiblock_grid grid_t,
        typename data_alias_t,
        dims::grid_array_dimension major_dim_t=dims::singleton_dim,
        const array_centering centering=cell_centered,
        array_containers::is_array_container container_t = default_container_t<typename detail::get_fundamental_type<data_alias_t>::type>
        >
    struct grid_array
    {
        constexpr static array_centering centering_type() { return centering; }
        static constexpr std::size_t total_idx_rank(void) {return minor_dim_t::rank()+dims::dynamic_dims<4>::rank()+major_dim_t::rank();}
        
        typedef grid_t grid_type;
        typedef typename detail::get_fundamental_type<data_alias_t>::type fundamental_type;
        typedef fundamental_type value_type;
        typedef typename detail::get_dim_type<data_alias_t>::type minor_dim_t;
        typedef minor_dim_t array_minor_dim_t;
        typedef major_dim_t array_major_dim_t;
        typedef typename detail::get_index_type<centering_type()>::integral_type index_integral_type;
        typedef typename detail::get_index_type<centering_type()>::array_type    index_type;
        typedef typename grid_t::coord_point_type coord_point_type;
        typedef data_alias_t alias_type;
        typedef data_alias_t unwrapped_minor_type;
        
        typedef typename detail::get_variable_mem_map<data_alias_t>::type variable_map_type;
        typedef typename detail::get_ijklb_map_type<centering_type(), grid_type::dim()>::type grid_map_type;
        typedef composite_map_t<variable_map_type, grid_map_type> mem_map_type;
        
        const grid_t* grid;
        minor_dim_t minor_dims;
        dims::dynamic_dims<detail::get_centering_grid_idx_rank<centering_type()>::value> grid_dims;
        major_dim_t major_dims;
        container_t data;
        std::size_t offset;
        std::size_t idx_coeffs [total_idx_rank()];
        mem_map_type mem_map;
        
        const auto& var_map() const
        {
            return std::get<0>(mem_map.maps);
        }
        
        const auto& block_map() const
        {
            return std::get<1>(mem_map.maps);
        }
        
        auto& var_map()
        {
            return std::get<0>(mem_map.maps);
        }
        
        auto& block_map()
        {
            return std::get<1>(mem_map.maps);
        }
        
        grid_array(){}
        grid_array(
            const grid_t& grid_in,
            const data_alias_t& fill_elem,
            const major_dim_t& major_dims_in=dims::singleton_dim())
        {
            grid = &grid_in;
            minor_dims = minor_dim_t();
            major_dims = major_dims_in;
            grid_dims = create_grid_dims<decltype(grid_in), centering>(grid_in);
            data.resize(minor_dims.total_size()*grid_dims.total_size()*major_dims.total_size());
            std::size_t n = total_idx_rank();
            for (auto i: range(0, n)) idx_coeffs[i] = 1;
            for (auto i: range(0, n))
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
            auto& grid_map = block_map();
            detail::insert_grid_dims<centering_type()>(grid_map, grid_in);
            this->mem_map.compute_coeffs();
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
        
        template <const int lev, typename idx_t>
        requires (std::convertible_to<idx_t, int>)
         _finline_ std::size_t offst_r(const idx_t& idx) const
        {
            return idx_coeffs[lev]*idx;
        }
        
        template <const int lev, typename idx_t, typename... idxs_t>
        requires (std::convertible_to<idx_t, int>)
        _finline_ std::size_t offst_r(const idx_t& idx, const idxs_t&... idxs) const
        {
            return idx_coeffs[lev]*idx + offst_r<lev+1>(idxs...);
        }
        
        //Note: This has potential to become a bottleneck. Consider a strategy of precomputing the zero-offset, then
        //indexing everything by simply multiplying individual indices by coefficients.
        //Also, try incrementing the computed offset indx-for-index instead of constructing properly-ranked arrays
        template <typename... idxs_t>
        _finline_ fundamental_type& operator() (const idxs_t&... idxs)
        {
            // static_assert(detail::index_rank_size_t<idxs_t...>::value == total_idx_rank(), "wrong number of indices passed to indexing operator");
            // return data[offset + offst_r<0>(idxs...)];
            return data[mem_map.offset(idxs...)];
        }
        
        //Note: This has potential to become a bottleneck.
        template <typename... idxs_t>
        _finline_ const fundamental_type& operator() (const idxs_t&... idxs) const
        {
            // static_assert(detail::index_rank_size_t<idxs_t...>::value == total_idx_rank(), "wrong number of indices passed to indexing operator");
            // return data[offset + offst_r<0>(idxs...)];
            return data[mem_map.offset(idxs...)];
        }
        
        void set_elem(const index_type& idx, const alias_type& alias)
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int i = 0; i < alias.size(); ++i) (*this)(i, idx) = alias[i];
            }
            if constexpr (!ctrs::basic_array<alias_type>)
            {
                (*this)(idx) = alias;
            }
        }
        
        alias_type get_elem(const index_type& idx) const
        {
            if constexpr (ctrs::basic_array<alias_type>)
            {
                alias_type output;
                for (int i = 0; i < output.size(); ++i) output[i] = (*this)(i, idx);
                return output;
            }
            if constexpr (!ctrs::basic_array<alias_type>)
            {
                return (*this)(idx);
            }
        }
        
        // _finline_ fundamental_type& operator () (const int i, const int j, const int k, const int lb)
        // {
        //     // std::size_t off =
        //     //     idx_coeffs[0]*(i-grid->get_num_cells(0)) +
        //     //     idx_coeffs[1]*(j-grid->get_num_cells(1)) +
        //     //     idx_coeffs[2]*(k-grid->get_num_cells(2)) +
        //     //     idx_coeffs[3]*(lb);
        //     std::size_t off = offset + 
        //         idx_coeffs[0]*(i-grid->get_num_cells(0)) +
        //         idx_coeffs[1]*(j-grid->get_num_cells(1)) +
        //         idx_coeffs[2]*(k-grid->get_num_cells(2)) +
        //         idx_coeffs[3]*(lb);
        //     return data[off];
        // }
        
        // _finline_ const fundamental_type& operator () (const int i, const int j, const int k, const int lb) const
        // {
        //     // std::size_t off =
        //     //     idx_coeffs[0]*(i-grid->get_num_cells(0)) +
        //     //     idx_coeffs[1]*(j-grid->get_num_cells(1)) +
        //     //     idx_coeffs[2]*(k-grid->get_num_cells(2)) +
        //     //     idx_coeffs[3]*(lb);
        //     std::size_t off = offset + 
        //         idx_coeffs[0]*(i) +
        //         idx_coeffs[1]*(j) +
        //         idx_coeffs[2]*(k) +
        //         idx_coeffs[3]*(lb);
        //     return data[off];
        // }
        
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
                i4*idx_coeffs[minor_dim_t::rank()+3]//+
                // i5*idx_coeffs[minor_dim_t::rank()+4]
            ];
        }
        
        template <typename i0_t, typename i1_t, typename i2_t, typename i3_t, typename i4_t, typename i5_t>
        requires std::convertible_to<i0_t, int> &&
        std::convertible_to<i1_t, int> &&
        std::convertible_to<i2_t, int> &&
        std::convertible_to<i3_t, int> &&
        std::convertible_to<i4_t, int> &&
        std::convertible_to<i5_t, int>
        _finline_ const fundamental_type& unwrap_idx(const i0_t& i0, const i1_t& i1, const i2_t& i2, const i3_t& i3, const i4_t& i4, const i5_t& i5) const
        {
            return data[
                offset +
                i0*idx_coeffs[0] + 
                i1*idx_coeffs[minor_dim_t::rank()+0]+
                i2*idx_coeffs[minor_dim_t::rank()+1]+
                i3*idx_coeffs[minor_dim_t::rank()+2]+
                i4*idx_coeffs[minor_dim_t::rank()+3]//+
                // i5*idx_coeffs[minor_dim_t::rank()+4]
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
            //Note: attempted to hand-optimize this, turns out the compiler will do just fine on its own.
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
        
        minor_dim_t get_minor_dims() const { return minor_dims; }
        major_dim_t get_major_dims() const { return major_dims; }
        major_dim_t get_grid_dims () const { return grid_dims;  }
        auto get_total_dims()        const { return minor_dims*grid_dims*major_dims; }
        const grid_t& get_grid() const {return *grid;}
        
        auto get_grid_index_bounding_box(const grid::exchange_inclusion_e& exchange_policy) const
        {
            const auto& grid = get_grid();
            int iexchg = 0;
            if (exchange_policy == include_exchanges) iexchg = 1;
            bound_box_t<int, grid_map_type::rank()> output;
            const auto& ijk_map = this->block_map();
            
            output.min(index_type::i_idx) = std::get<index_type::i_idx>(ijk_map.dims).i0 + iexchg*grid.get_num_exchange(0);
            output.max(index_type::i_idx) = std::get<index_type::i_idx>(ijk_map.dims).i1 - iexchg*grid.get_num_exchange(0);
            
            output.min(index_type::j_idx) = std::get<index_type::j_idx>(ijk_map.dims).i0 + iexchg*grid.get_num_exchange(1);
            output.max(index_type::j_idx) = std::get<index_type::j_idx>(ijk_map.dims).i1 - iexchg*grid.get_num_exchange(1);
    
            output.min(index_type::k_idx) = std::get<index_type::k_idx>(ijk_map.dims).i0 + iexchg*grid.get_num_exchange(2);
            output.max(index_type::k_idx) = std::get<index_type::k_idx>(ijk_map.dims).i1 - iexchg*grid.get_num_exchange(2);
            
            output.min(index_type::lb_idx) = std::get<index_type::lb_idx>(ijk_map.dims).i0;
            output.max(index_type::lb_idx) = std::get<index_type::lb_idx>(ijk_map.dims).i1;

            detail::set_bbox_dir_val<centering_type()>(output, grid_type::dim());
            
            return output;
        }
    };
    
    template <
        multiblock_grid grid_t,
        typename data_alias_t,
        dims::grid_array_dimension major_dim_t = dims::singleton_dim,
        array_containers::is_array_container container_t = default_container_t<typename detail::get_fundamental_type<data_alias_t>::type>
        > using face_array = grid_array<grid_t, data_alias_t, major_dim_t, face_centered, container_t>;
    template <
        multiblock_grid grid_t,
        typename data_alias_t,
        dims::grid_array_dimension major_dim_t = dims::singleton_dim,
        array_containers::is_array_container container_t = default_container_t<typename detail::get_fundamental_type<data_alias_t>::type>
        > using cell_array = grid_array<grid_t, data_alias_t, major_dim_t, cell_centered, container_t>;
}