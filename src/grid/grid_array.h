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
#include "grid/partition.h"
#include "core/static_math.h"
#include "grid/grid_index_types.h"
#include "grid/grid.h"
#include "grid/array_image.h"
#include "core/mem_map.h"
#include "core/vec_image.h"
#include "core/permute.h"

#include "dispatch/device_type.h"
#include "dispatch/device_vector.h"
#include "dispatch/execute.h"
#include "dispatch/ranges/linear_range.h"
#include "dispatch/support_of.h"

#include "algs/transform_inplace.h"


namespace spade::grid
{
    namespace detail
    {        
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
        
        template <const array_centering centering, typename grid_map_t, typename grid_t, typename iarray_t>
        requires (centering == cell_centered)
        void insert_grid_dims(grid_map_t& map, const grid_t& grid, const iarray_t& exch)
        {
            //i
            std::get<cell_idx_t::i_idx>(map.views).i0 = -exch[0];
            std::get<cell_idx_t::i_idx>(map.views).i1 =  grid.get_num_cells(0) + exch[0];
            
            //j
            std::get<cell_idx_t::j_idx>(map.views).i0 = -exch[1];
            std::get<cell_idx_t::j_idx>(map.views).i1 =  grid.get_num_cells(1) + exch[1];
            
            //k
            std::get<cell_idx_t::k_idx>(map.views).i0 = -exch[2];
            std::get<cell_idx_t::k_idx>(map.views).i1 =  grid.get_num_cells(2) + exch[2];
            
            //lb
            std::get<cell_idx_t::lb_idx>(map.views).i0 = 0;
            std::get<cell_idx_t::lb_idx>(map.views).i1 = grid.get_num_local_blocks();
        }
        
        template <const array_centering centering, typename grid_map_t, typename grid_t, typename iarray_t>
        requires (centering == face_centered)
        void insert_grid_dims(grid_map_t& map, const grid_t& grid, const iarray_t& exch)
        {
            //i
            std::get<face_idx_t::i_idx>(map.views).i0 = -exch[0];
            std::get<face_idx_t::i_idx>(map.views).i1 =  grid.get_num_cells(0) + exch[0] + 1;
            
            //j
            std::get<face_idx_t::j_idx>(map.views).i0 = -exch[1];
            std::get<face_idx_t::j_idx>(map.views).i1 =  grid.get_num_cells(1) + exch[1] + 1;
            
            //k
            std::get<face_idx_t::k_idx>(map.views).i0 = -exch[2];
            std::get<face_idx_t::k_idx>(map.views).i1 =  grid.get_num_cells(2) + exch[2] + 1;
            
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
    
    template <
        multiblock_grid grid_t,
        typename        data_alias_t,
        typename        device_t,
        typename        mmap_t = mem_map::tlinear_t,
        const           array_centering centering = cell_centered
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
        
        template <typename c_data_t>
        using container_type
            = std::conditional<
                device::is_cpu<device_t>,
                std::vector<c_data_t>,
                device::device_vector<c_data_t>
                >::type;
        
        using variable_map_type   = detail::get_variable_mem_map<data_alias_t>::type;
        using grid_map_type       = detail::get_ijklb_map_type<centering_type(), grid_type::dim()>::type;
        using device_type         = device_t;
        using mem_map_type        = decltype(mem_map::make_grid_map(mmap_t(), alias_type(), device_type(), {0, 0}, {0, 0}, {0, 0}, {0, 0}));
        using var_idx_t           = ctrs::array<int, variable_map_type::rank()>;
        
        using image_type           = array_image_t<alias_type, mem_map_type,       value_type*, centering_type(), grid_type>;
        using const_image_type     = array_image_t<alias_type, mem_map_type, const value_type*, centering_type(), grid_type>;
        using iarray_t             = ctrs::array<int, 3>;
        
        iarray_t                         num_exch;
        const grid_t*                    grid;
        container_type<fundamental_type> data;
        mem_map_type                     mmap;
        std::size_t                      var_offset;
        
        constexpr static int dim() { return grid_t::dim(); }
        
        static_assert(!(device::is_gpu<device_t> && !_sp_cuda), "attempted to declare GPU array without GPU support");
        
        device_t device() const { return device_t(); }
        
        fundamental_type* get_base() { return &data[0]; }
        
        void clear()
        {
            data.clear();
        }
        
        std::size_t get_base_offset(const index_type& ii) const
        {
            return mmap.compute_offset(var_idx_t(0), ii);
        }
        
        int      get_num_exchange(int i) const { return num_exch[i]; }
        iarray_t get_num_exchange()      const { return num_exch; }
        
        auto get_device() const {return device_t();}

        // std::size_t cell_element_size() const {return mem_map::map_size(var_map());}

        grid_array(){}
        grid_array(
            const grid_t& grid_in,
            const data_alias_t& fill_elem,
            const typename grid_t::array_desig_type& num_exch_in,
            const device_t& dev_in,
            const mmap_t& mmap_tag = mem_map::linear
            ) : grid{&grid_in}, num_exch{num_exch_in}
        {
            mmap = mem_map::make_grid_map(
                mmap_tag, alias_type(), device_type(),
                ctrs::make_array(-num_exch[0], grid->get_num_cells(0) + num_exch[0]),
                ctrs::make_array(-num_exch[1], grid->get_num_cells(1) + num_exch[1]),
                ctrs::make_array(-num_exch[2], grid->get_num_cells(2) + num_exch[2]),
                ctrs::make_array(0, int(grid->get_num_local_blocks())));
            auto total_size = mmap.volume();
            data.resize(total_size);
        }
        
        grid_array(const grid_array& rhs) = default;
        
        grid_array(grid_array&& rhs)
        {
            grid     = rhs.grid;
            num_exch = rhs.num_exch;
            data     = std::move(rhs.data);
            mmap     = rhs.mmap;
        }
        
        grid_array& operator = (const grid_array& rhs) = default;
        
        grid_array& operator = (grid_array&& rhs)
        {
            grid     = rhs.grid;
            num_exch = rhs.num_exch;
            data     = std::move(rhs.data);
            mmap     = rhs.mmap;
            return *this;
        }
        
        const const_image_type image() const { return {&data[0], data.size(), mmap, var_offset}; }
        image_type             image()       { return {&data[0], data.size(), mmap, var_offset}; }

        //TODO: clean this up a little bit with more lambdas!
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator -= (const rhs_t& rhs)
        {
            auto my_img = this->image();
            const auto rhs_img = rhs.image();
            const auto rg = dispatch::support_of(*this, grid::exclude_exchanges);
            
            auto loop = [=] _sp_hybrid (const index_type& idx) mutable
            {
                const auto rhs_elem = rhs_img.get_elem(idx);
                my_img.decr_elem(idx, rhs_elem);
            };
            
            dispatch::execute(rg, loop);
            return *this;
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator += (const rhs_t& rhs)
        {
            auto my_img = this->image();
            const auto rhs_img = rhs.image();
            const auto rg = dispatch::support_of(*this, grid::exclude_exchanges);
            
            auto loop = [=] _sp_hybrid (const index_type& idx) mutable
            {
                const auto rhs_elem = rhs_img.get_elem(idx);
                my_img.incr_elem(idx, rhs_elem);
            };
            
            dispatch::execute(rg, loop);
            return *this;
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator *= (const rhs_t& rhs)
        {
            auto my_img = this->image();
            const auto rhs_img = rhs.image();
            const auto rg = dispatch::support_of(*this, grid::exclude_exchanges);
            
            auto loop = [=] _sp_hybrid (const index_type& idx) mutable
            {
                const auto rhs_elem = rhs_img.get_elem(idx);
                my_img.mult_elem(idx, rhs_elem);
            };
            
            dispatch::execute(rg, loop);
            return *this;
        }
        
        template <multiblock_array rhs_t>
        requires elementwise_compatible<grid_array, rhs_t>
        grid_array& operator /= (const rhs_t& rhs)
        {
            auto my_img = this->image();
            const auto rhs_img = rhs.image();
            const auto rg = dispatch::support_of(*this, grid::exclude_exchanges);
            
            auto loop = [=] _sp_hybrid (const index_type& idx) mutable
            {
                const auto rhs_elem = rhs_img.get_elem(idx);
                my_img.divi_elem(idx, rhs_elem);
            };
            
            dispatch::execute(rg, loop);
            return *this;
        }
        
        template <typename numeric_t> grid_array& operator *= (const numeric_t& rhs)
        requires std::floating_point<numeric_t>
        {
            std::size_t i0 = 0;
            std::size_t i1 = this->data.size();
            const dispatch::ranges::linear_range_t lrange(i0, i1, this->device());
            auto vimg         = utils::make_vec_image(this->data);
            auto fnc = [=] _sp_hybrid (const std::size_t& i) mutable { vimg[i] *= rhs; };
            dispatch::execute(lrange, fnc);
            return *this;
        }
        
        template <typename numeric_t> grid_array& operator /= (const numeric_t& rhs)
        requires std::floating_point<numeric_t>
        {
            std::size_t i0 = 0;
            std::size_t i1 = this->data.size();
            const dispatch::ranges::linear_range_t lrange(i0, i1, this->device());
            auto vimg         = utils::make_vec_image(this->data);
            auto fnc = [=] _sp_hybrid (const std::size_t& i) mutable { vimg[i] /= rhs; };
            dispatch::execute(lrange, fnc);
            return *this;
        }
        
        grid_array& operator = (const fundamental_type& rhs)
        {
            std::size_t i0 = 0;
            std::size_t i1 = this->data.size();
            const dispatch::ranges::linear_range_t lrange(i0, i1, this->device());
            auto vimg         = utils::make_vec_image(this->data);
            auto fnc = [=] _sp_hybrid (const std::size_t& i) mutable { vimg[i] = rhs; };
            dispatch::execute(lrange, fnc);
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
}