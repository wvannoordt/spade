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
#include "core/partition.h"
#include "core/static_math.h"
#include "core/grid_index_types.h"
#include "core/grid.h"
#include "core/mem_map.h"

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
            t.size(); //todo: figure out whay the heck is going on with this
        };
                
        template <typename data_t> struct get_variable_mem_map
        {
            typedef mem_map::empty_view_t type;
        };
        
        template <is_static_1D_array data_t> struct get_variable_mem_map<data_t>
        {
            typedef mem_map::recti_view_t<mem_map::static_dim_t<0,data_t::size()>> type;
        };
        
        
        //getting the memory map type for the grid indices
        template <const array_centering center, const std::size_t grid_dim> struct get_ijklb_map_type{};
        
        template <const std::size_t grid_dim> struct get_ijklb_map_type<cell_centered, grid_dim>
        {
            typedef mem_map::recti_view_t<
                mem_map::dynamic_dim_t<int>,
                mem_map::dynamic_dim_t<int>,
                mem_map::dynamic_dim_t<int>,
                mem_map::dynamic_dim_t<int>> type;
        };
        
        template <const std::size_t grid_dim> struct get_ijklb_map_type<node_centered, grid_dim>
        {
            typedef mem_map::recti_view_t<
                mem_map::dynamic_dim_t<int>,
                mem_map::dynamic_dim_t<int>,
                mem_map::dynamic_dim_t<int>,
                mem_map::dynamic_dim_t<int>> type;
        };
        
        template <const std::size_t grid_dim> struct get_ijklb_map_type<face_centered, grid_dim>
        {
            //Can I do a nicer job here?
            typedef mem_map::recti_view_t<
                typename std::conditional<0==face_idx_t::dir_idx, mem_map::static_dim_t<0, grid_dim>, mem_map::dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, mem_map::static_dim_t<0, grid_dim>, mem_map::dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, mem_map::static_dim_t<0, grid_dim>, mem_map::dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, mem_map::static_dim_t<0, grid_dim>, mem_map::dynamic_dim_t<int>>::type,
                typename std::conditional<0==face_idx_t::dir_idx, mem_map::static_dim_t<0, grid_dim>, mem_map::dynamic_dim_t<int>>::type
                > type;
        };
        
        template <const array_centering centering> struct get_centering_grid_idx_rank
        {
            static constexpr int value = 4;
        };
        
        template <> struct get_centering_grid_idx_rank<cell_centered> { static constexpr int value = 4; };
        template <> struct get_centering_grid_idx_rank<face_centered> { static constexpr int value = 5; };
        template <> struct get_centering_grid_idx_rank<node_centered> { static constexpr int value = 4; };
        
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
            std::get<cell_idx_t::i_idx>(map.views).i0 = -grid.get_num_exchange(0);
            std::get<cell_idx_t::i_idx>(map.views).i1 =  grid.get_num_cells(0) + grid.get_num_exchange(0);
            
            //j
            std::get<cell_idx_t::j_idx>(map.views).i0 = -grid.get_num_exchange(1);
            std::get<cell_idx_t::j_idx>(map.views).i1 =  grid.get_num_cells(1) + grid.get_num_exchange(1);
            
            //k
            std::get<cell_idx_t::k_idx>(map.views).i0 = -grid.get_num_exchange(2);
            std::get<cell_idx_t::k_idx>(map.views).i1 =  grid.get_num_cells(2) + grid.get_num_exchange(2);
            
            //lb
            std::get<cell_idx_t::lb_idx>(map.views).i0 = 0;
            std::get<cell_idx_t::lb_idx>(map.views).i1 = grid.get_num_local_blocks();
        }
        
        template <const array_centering centering, typename grid_map_t, typename grid_t>
        requires (centering == face_centered)
        void insert_grid_dims(grid_map_t& map, const grid_t& grid)
        {
            //i
            std::get<face_idx_t::i_idx>(map.views).i0 = -grid.get_num_exchange(0);
            std::get<face_idx_t::i_idx>(map.views).i1 =  grid.get_num_cells(0) + grid.get_num_exchange(0) + 1;
            
            //j
            std::get<face_idx_t::j_idx>(map.views).i0 = -grid.get_num_exchange(1);
            std::get<face_idx_t::j_idx>(map.views).i1 =  grid.get_num_cells(1) + grid.get_num_exchange(1) + 1;
            
            //k
            std::get<face_idx_t::k_idx>(map.views).i0 = -grid.get_num_exchange(2);
            std::get<face_idx_t::k_idx>(map.views).i1 =  grid.get_num_cells(2) + grid.get_num_exchange(2) + 1;
            
            //lb
            std::get<face_idx_t::lb_idx>(map.views).i0 = 0;
            std::get<face_idx_t::lb_idx>(map.views).i1 = grid.get_num_local_blocks();
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
        const array_centering centering=cell_centered,
        array_containers::is_array_container container_t = default_container_t<typename detail::get_fundamental_type<data_alias_t>::type>
        >
    struct grid_array
    {
        constexpr static array_centering centering_type() { return centering; }
        
        using grid_type           = grid_t;
        using fundamental_type    = detail::get_fundamental_type<data_alias_t>::type;
        using value_type          = fundamental_type;
        using index_integral_type = get_index_type<centering_type()>::integral_type;
        using index_type          = get_index_type<centering_type()>::array_type;
        using coord_point_type    = grid_t::coord_point_type;
        using alias_type          = data_alias_t;
        
        using variable_map_type   = detail::get_variable_mem_map<data_alias_t>::type;
        using grid_map_type       = detail::get_ijklb_map_type<centering_type(), grid_type::dim()>::type;
        using mem_map_type        = mem_map::mem_map_t<mem_map::recti_view_t<variable_map_type, grid_map_type>>;
        
        const grid_t* grid;
        container_t data;
        std::size_t offset;
        mem_map_type mem_view;
        
        const auto& var_map() const
        {
            return std::get<0>(mem_view.mmap.views);
        }
        
        const auto& block_map() const
        {
            return std::get<1>(mem_view.mmap.views);
        }
        
        auto& var_map()
        {
            return std::get<0>(mem_view.mmap.views);
        }
        
        auto& block_map()
        {
            return std::get<1>(mem_view.mmap.views);
        }

        std::size_t cell_element_size() const {return mem_map::map_size(var_map());}

        grid_array(){}
        grid_array(
            const grid_t& grid_in,
            const data_alias_t& fill_elem
            )
        {
            grid = &grid_in;
            auto& grid_map = block_map();
            detail::insert_grid_dims<centering_type()>(grid_map, grid_in);
            this->mem_view.compute_coeffs();
            this->mem_view.compute_offset_base();
            auto total_size = mem_map::map_size(var_map())*mem_map::map_size(block_map());
            data.resize(total_size);
        }

        template <typename... idxs_t>
        _finline_ fundamental_type& operator() (const idxs_t&... idxs)
        {
            return data[mem_view.compute_offset(idxs...)];
        }

        //Note: This has potential to become a bottleneck.
        template <typename... idxs_t>
        _finline_ const fundamental_type& operator() (const idxs_t&... idxs) const
        {
            return data[mem_view.compute_offset(idxs...)];
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

        const grid_t& get_grid() const {return *grid;}
        
        auto get_grid_index_bounding_box(const grid::exchange_inclusion_e& exchange_policy) const
        {
            const auto& grid = get_grid();
            int iexchg = 1;
            if (exchange_policy == include_exchanges) iexchg = 0;
            bound_box_t<int, grid_map_type::rank()> output;
            const auto& ijk_map = this->block_map();
            
            output.min(index_type::i_idx) = std::get<index_type::i_idx>(ijk_map.views).i0 + iexchg*grid.get_num_exchange(0);
            output.max(index_type::i_idx) = std::get<index_type::i_idx>(ijk_map.views).i1 - iexchg*grid.get_num_exchange(0);
            
            output.min(index_type::j_idx) = std::get<index_type::j_idx>(ijk_map.views).i0 + iexchg*grid.get_num_exchange(1);
            output.max(index_type::j_idx) = std::get<index_type::j_idx>(ijk_map.views).i1 - iexchg*grid.get_num_exchange(1);
    
            output.min(index_type::k_idx) = std::get<index_type::k_idx>(ijk_map.views).i0 + iexchg*grid.get_num_exchange(2);
            output.max(index_type::k_idx) = std::get<index_type::k_idx>(ijk_map.views).i1 - iexchg*grid.get_num_exchange(2);
            
            output.min(index_type::lb_idx) = std::get<index_type::lb_idx>(ijk_map.views).i0;
            output.max(index_type::lb_idx) = std::get<index_type::lb_idx>(ijk_map.views).i1;
            
            detail::set_bbox_dir_val<centering_type()>(output, grid_type::dim());
            
            return output;
        }
    };
    
    template <
        multiblock_grid grid_t,
        typename data_alias_t,
        array_containers::is_array_container container_t = default_container_t<typename detail::get_fundamental_type<data_alias_t>::type>
        > using face_array = grid_array<grid_t, data_alias_t, face_centered, container_t>;
    template <
        multiblock_grid grid_t,
        typename data_alias_t,
        array_containers::is_array_container container_t = default_container_t<typename detail::get_fundamental_type<data_alias_t>::type>
        > using cell_array = grid_array<grid_t, data_alias_t, cell_centered, container_t>;
}