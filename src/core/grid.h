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

namespace spade::grid
{    
    template <class T> concept multiblock_grid = requires(T t, const cell_idx_t& i_c, const face_idx_t& i_f, const node_idx_t& i_n)
    {
        // todo: write this
        { t.get_coords(i_c) } -> ctrs::basic_array;
        { t.get_coords(i_f) } -> ctrs::basic_array;
        { t.get_coords(i_n) } -> ctrs::basic_array;
        { t.get_comp_coords(i_c) } -> ctrs::basic_array;
        { t.get_comp_coords(i_f) } -> ctrs::basic_array;
        { t.get_comp_coords(i_n) } -> ctrs::basic_array;
    };
    
    
    enum array_center_e
    {
        cell_centered=0,
        node_centered=1
    };
    
    template <class T, const array_center_e ct> concept has_centering_type = (T::centering_type() == ct);
    
    template <class T> concept multiblock_array = requires(T t, int a, int i, int j, int k, int lb, int b)
    {
        { t.get_major_dims() } -> dims::grid_array_dimension;
        { t.get_minor_dims() } -> dims::grid_array_dimension;
        { t.get_grid_dims()  } -> dims::grid_array_dimension;
        { t.get_grid()       } -> multiblock_grid;
        t.unwrap_idx(a, i, j, k, lb, b);
    };
    
    template <typename T0, typename T1> concept elementwise_compatible = multiblock_array<T0> && multiblock_array<T1> &&
    requires(const T0& t0, const T1& t1)
    {
        //TODO: write this
        t0;
    };
    
    enum exchange_inclusion_e
    {
        exclude_exchanges=0,
        include_exchanges=1
    };
    
    template <multiblock_grid grid_t> static auto create_grid_dims(const grid_t& grid, const array_center_e& center_type)
    {
        switch (center_type)
        {
            case cell_centered:
            {
                return dims::dynamic_dims<4>(
                    grid.get_num_cells(0)+grid.get_num_exchange(0)*2,
                    grid.get_num_cells(1)+grid.get_num_exchange(1)*2,
                    grid.get_num_cells(2)+grid.get_num_exchange(2)*2,
                    grid.get_num_local_blocks());
            }
            case node_centered:
            {
                return dims::dynamic_dims<4>(
                    1+grid.get_num_cells(0)+grid.get_num_exchange(0)*2,
                    1+grid.get_num_cells(1)+grid.get_num_exchange(1)*2,
                    1+grid.get_num_cells(2)+grid.get_num_exchange(2)*2,
                    grid.get_num_local_blocks());
            }
        }
        return dims::dynamic_dims<4>(0,0,0,0);
    }
    
    struct neighbor_relationship_t
    {
        ctrs::array<int, 3> edge_vec;
        int rank_end, rank_start;
        int lb_glob_end, lb_glob_start;
    };
    
    template <coords::coordinate_system coord_t, parallel::parallel_group par_group_t, typename dimension_t=static_math::int_const_t<3>>
    class cartesian_grid_t
    {
        public:
            typedef coord_t::coord_type dtype;
            typedef coord_t::coord_type coord_type;
            typedef coord_t coord_sys_type;
            typedef dimension_t dim_t;
            cartesian_grid_t(
                const ctrs::array<int,   dimension_t::value>& num_blocks_in,
                const ctrs::array<int,   dimension_t::value>& cells_in_block_in,
                const ctrs::array<int,   dimension_t::value>& exchange_cells_in,
                const bound_box_t<dtype, dimension_t::value>& bounds_in,
                const coord_t& coord_system_in,
                par_group_t& group_in,
                const dimension_t& dim_v = static_math::int_const_t<3>())
            {
                //Initialize class members
                this->grid_group = &group_in;
                coord_system = coord_system_in;
                dx = 1.0;
                num_blocks = 1;
                cells_in_block = 1;
                total_blocks = 1;
                is3d = (dim()==3);
                bounds.min(2) = 0.0;
                bounds.max(2) = 1.0;
                
                //copy constructor arguments, compute domain bounding blox
                for (std::size_t i = 0; i < dim(); i++)
                {
                    bounds.max(i) = bounds_in.max(i);
                    bounds.min(i) = bounds_in.min(i);
                    num_blocks[i] = num_blocks_in[i];
                    cells_in_block[i] = cells_in_block_in[i];
                    dx[i] = bounds.size(i) / (num_blocks[i]*cells_in_block[i]);
                    exchange_cells[i] = exchange_cells_in[i];
                    total_blocks *= num_blocks[i];
                }
                
                
                //partition the domain
                grid_partition = partition::block_partition_t(num_blocks, &group_in);
                
                //compute individual block bounding boxes
                block_boxes.resize(grid_partition.get_num_local_blocks());
                for (auto lb: range(0,grid_partition.get_num_local_blocks()))
                {
                    auto& box = block_boxes[lb];
                    ctrs::array<std::size_t, 3> glob_block_idx = ctrs::expand_index(grid_partition.get_global_block(lb), num_blocks);
                    static_for<0,dim()>([&](auto i)
                    {
                        box.min(i.value) = bounds.min(i.value) + (glob_block_idx[i.value]+0)*bounds.size(i.value)/num_blocks[i.value];
                        box.max(i.value) = bounds.min(i.value) + (glob_block_idx[i.value]+1)*bounds.size(i.value)/num_blocks[i.value];
                    });
                }
                
                
                block_is_domain_boundary.resize(this->get_num_global_blocks());
                for (auto lb: range(0,this->get_num_global_blocks()))
                {
                    const auto& lbi = lb;
                    auto& data = block_is_domain_boundary[lbi];
                    const auto lb_idx = ctrs::expand_index(lbi, num_blocks);
                    for (auto& val: data) val = false;
                    for (auto d: range(0,dim()))
                    {
                        data[2*d+0] = (lb_idx[d]==0);
                        data[2*d+1] = (lb_idx[d]==(num_blocks[d]-1));
                    }
                }
                
                send_size_elems.resize(group_in.size(), 0);
                recv_size_elems.resize(group_in.size(), 0);
                
                send_bufs.resize(group_in.size());
                recv_bufs.resize(group_in.size());
                
                requests.resize(group_in.size());
                statuses.resize(group_in.size());
                
                //compute neighbor relationships
                for (auto lb: range(0, this->get_num_global_blocks()))
                {
                    std::size_t lb_glob = lb;
                    ctrs::array<int, 3> delta_lb = 0;
                    int i3d = this->is_3d()?1:0;
                    auto delta_lb_range = range(-1, 2)*range(-1, 2)*range(-i3d, 1+i3d);
                    for (auto dlb: delta_lb_range)
                    {
                        ctrs::copy_array(dlb, delta_lb);
                        ctrs::array<int, 3> lb_nd;
                        ctrs::copy_array(expand_index(lb_glob, this->num_blocks), lb_nd);
                        lb_nd += delta_lb;
                        lb_nd += this->num_blocks;
                        lb_nd %= this->num_blocks;
                        std::size_t lb_glob_neigh = ctrs::collapse_index(lb_nd, this->num_blocks);                        
                        int rank_here  = grid_partition.get_global_rank(lb_glob);
                        int rank_neigh = grid_partition.get_global_rank(lb_glob_neigh);
                        
                        neighbor_relationship_t neighbor_relationship;
                        neighbor_relationship.edge_vec      = delta_lb;
                        neighbor_relationship.rank_end      = rank_neigh;
                        neighbor_relationship.lb_glob_end   = lb_glob_neigh;
                        neighbor_relationship.rank_start    = rank_here;
                        neighbor_relationship.lb_glob_start = lb_glob;
                        if (neighbor_relationship.edge_vec[0] != 0 || neighbor_relationship.edge_vec[1] != 0 || neighbor_relationship.edge_vec[2] != 0)
                        {
                            if ((rank_here==group_in.rank()) || (rank_neigh==group_in.rank()))
                            {
                                neighbors.push_back(neighbor_relationship);
                                
                                //Here sends, neigh receives
                                if (group_in.rank()==rank_here)
                                {
                                    send_size_elems[rank_neigh] += this->get_send_index_bounds(delta_lb).volume();
                                }
                                if (group_in.rank()==rank_neigh)
                                {
                                    recv_size_elems[rank_here]  += this->get_recv_index_bounds(delta_lb).volume();
                                }
                            }
                        }
                    }
                }
            }
            
            constexpr static int dim(void) {return dim_t::value;}
            
            template <multiblock_grid_idx_t idx_t>_finline_ ctrs::array<dtype, 3> get_coords(const idx_t& i) const
            {
                return coord_system.map(this->get_comp_coords(i));
            }
            
            template <multiblock_grid_idx_t idx_t>_finline_ ctrs::array<dtype, 3> get_comp_coords(const idx_t& i) const
            {
                const auto  idx_r = get_index_coord(i);
                const int   lb    = (int)i[3];
                const auto& bnd   = block_boxes[lb];
                ctrs::array<dtype, 3> output(0.0, 0.0, 0.0);
                for (auto d: range(0,dim()))
                {
                    output[d] = bnd.min(d)+idx_r[d]*this->get_dx(d);
                }
                return output;
            }
            
            const par_group_t& group(void) const
            {
                return *grid_group;
            }
            
            const coord_t& coord_sys(void) const {return coord_system;}
            
            md_range_t<int,4> get_range(const array_center_e& centering_in, const exchange_inclusion_e& do_guards=exclude_exchanges) const
            {
                int iexchg = 0;
                int i3d = 0;
                if (dim()==3) i3d = 1;
                if (do_guards==include_exchanges) iexchg = 1;
                switch (centering_in)
                {
                    case cell_centered:
                    {
                        return md_range_t<int,4>(
                            -iexchg*exchange_cells[0],cells_in_block[0]+iexchg*exchange_cells[0],
                            -iexchg*exchange_cells[1],cells_in_block[1]+iexchg*exchange_cells[1],
                            -i3d*iexchg*exchange_cells[2],cells_in_block[2]+i3d*iexchg*exchange_cells[2],
                            0,grid_partition.get_num_local_blocks());
                    }
                    case node_centered:
                    {
                        return md_range_t<int,4>(
                            -iexchg*exchange_cells[0],1+cells_in_block[0]+iexchg*exchange_cells[0],
                            -iexchg*exchange_cells[1],1+cells_in_block[1]+iexchg*exchange_cells[1],
                            -i3d*iexchg*exchange_cells[2],(1-i3d)+cells_in_block[2]+i3d*iexchg*exchange_cells[2],
                            0,grid_partition.get_num_local_blocks());
                    }
                    default: return md_range_t<int,4>(0,0,0,0,0,0,0,0);
                }
            }
            
            std::size_t get_num_interior_cells(void) const {return get_num_cells(0)*get_num_cells(1)*get_num_cells(2)*get_num_local_blocks();}
            
            std::size_t get_num_blocks(const std::size_t& i)   const {return num_blocks[i];}
            ctrs::array<std::size_t, 3> get_num_blocks(void)   const {return num_blocks;}
            std::size_t  get_num_local_blocks(void) const {return grid_partition.get_num_local_blocks();}
            std::size_t get_num_global_blocks(void) const {return num_blocks[0]*num_blocks[1]*num_blocks[2];}
            int get_num_cells(const std::size_t& i)    const {return cells_in_block[i];}
            int get_num_exchange(const std::size_t& i) const {return exchange_cells[i];}
            bool is_domain_boundary(const std::size_t& lb_glob, const std::size_t& dir) const {return block_is_domain_boundary[lb_glob][dir];}
            bound_box_t<dtype,  3> get_bounds(void) const {return bounds;}
            ctrs::array<dtype,  3> get_dx(void) const {return dx;}
            dtype get_dx(const std::size_t& i) const {return dx[i];}
            bound_box_t<dtype, 3> get_block_box(const std::size_t& lb) const {return block_boxes[lb];}
            bool is_3d(void) const {return this->is3d;}
            std::size_t collapse_block_num(const std::size_t& lb_i, const std::size_t& lb_j, const std::size_t& lb_k) const
            {
                return lb_i + lb_j*this->get_num_blocks(0) + lb_k*this->get_num_blocks(0)*this->get_num_blocks(1);
            }
            //lb + edge_vec = lb_neigh (edge_vec points towards neighbor)
            bound_box_t<int, 3> get_recv_index_bounds(const ctrs::array<int, 3>& edge_vec) const
            {
                bound_box_t<int, 3> output;
                output.min(0) = 0; output.max(0) = 1;
                output.min(1) = 0; output.max(1) = 1;
                output.min(2) = 0; output.max(2) = 1;
                for (int i = 0; i < dim(); i++)
                {
                    ctrs::array<int, 3> recv_start(
                        (int)cells_in_block[i],
                        0,
                        -((int)exchange_cells[i]));
                    ctrs::array<int, 3> recv_end(
                        (int)(cells_in_block[i] + exchange_cells[i]),
                        (int)cells_in_block[i],
                        0);
                    output.min(i) = recv_start[1+edge_vec[i]];
                    output.max(i) = recv_end  [1+edge_vec[i]];
                }
                return output;
            }
            
            bound_box_t<int, 3> get_send_index_bounds(const ctrs::array<int, 3>& edge_vec) const
            {
                bound_box_t<int, 3> output;
                output.min(0) = 0; output.max(0) = 1;
                output.min(1) = 0; output.max(1) = 1;
                output.min(2) = 0; output.max(2) = 1;
                for (int i = 0; i < dim(); i++)
                {
                    ctrs::array<int, 3> send_start(
                        0,
                        0,
                        (int)cells_in_block[i]-(int)exchange_cells[i]);
                    ctrs::array<int, 3> send_end(
                        (int)exchange_cells[i],
                        (int)cells_in_block[i],
                        (int)cells_in_block[i]);
                    output.min(i) = send_start[1+edge_vec[i]];
                    output.max(i) = send_end  [1+edge_vec[i]];
                }
                return output;
            }
            
            template <grid::multiblock_array array_t> void exchange_array(array_t& array)
            {
                const std::size_t num_neighs = static_math::pow<3,dim()>::value;
                
                typedef typename array_t::value_type array_data_t;
                
                std::vector<std::size_t> num_send_bytes(grid_group->size(), 0);
                std::vector<std::size_t> num_recv_bytes(grid_group->size(), 0);
                
                std::size_t cell_elem_size = sizeof(array_data_t)*array.get_minor_dims().total_size();
                
                for (auto p:range(0, grid_group->size()))
                {
                    const std::size_t total_send_buf_size = cell_elem_size*send_size_elems[p];
                    const std::size_t total_recv_buf_size = cell_elem_size*recv_size_elems[p];
                    if (send_bufs[p].size() < total_send_buf_size) send_bufs[p].resize(total_send_buf_size);
                    if (recv_bufs[p].size() < total_recv_buf_size) recv_bufs[p].resize(total_recv_buf_size);
                }
                
                std::vector<std::size_t> offsets(grid_group->size(), 0);
                
                //pack
                for (auto& neigh_edge: neighbors)
                {
                    for (auto maj: range(0,array.get_major_dims().total_size()))
                    {
                        if (neigh_edge.rank_start==grid_group->rank())
                        {
                            std::vector<char>& send_buf_loc = send_bufs[neigh_edge.rank_end];
                            std::size_t& send_offset_loc = offsets[neigh_edge.rank_end];
                            std::size_t lb_loc = grid_partition.get_local_block(neigh_edge.lb_glob_start);
                            auto bounds = this->get_send_index_bounds(neigh_edge.edge_vec);
                            std::size_t copy_size = bounds.size(0)*cell_elem_size;                            
                            for (int k = bounds.min(2); k < bounds.max(2); ++k)
                            {
                                for (int j = bounds.min(1); j < bounds.max(1); ++j)
                                {
                                    std::copy(
                                        (char*)(&array.unwrap_idx(0,bounds.min(0),j,k,lb_loc,maj)),
                                        (char*)(&array.unwrap_idx(0,bounds.max(0),j,k,lb_loc,maj)),
                                        (char*)(&(send_buf_loc[send_offset_loc])));
                                    send_offset_loc += copy_size;
                                }
                            }
                        }
                    }
                }
                
                for (std::size_t p = 0; p < grid_group->size(); ++p)
                {
                    request_t r = grid_group->async_recv(&recv_bufs[p][0], recv_bufs[p].size(), p);
                    requests[p] = r;
                }
                for (std::size_t p = 0; p < grid_group->size(); ++p)
                {
                    grid_group->sync_send(&send_bufs[p][0], send_bufs[p].size(), p);
                }
                
                //might be a bottleneck
                grid_group->await_all(statuses.size(), requests.data(), statuses.data());
                
                //reset offsets
                for (auto& offset: offsets) offset = 0;
                
                //unpack
                for (auto& neigh_edge: neighbors)
                {
                    for (auto maj: range(0,array.get_major_dims().total_size()))
                    {
                        if (neigh_edge.rank_end==grid_group->rank())
                        {
                            std::vector<char>& recv_buf_loc = recv_bufs[neigh_edge.rank_start];
                            std::size_t& recv_offset_loc = offsets[neigh_edge.rank_start];
                            std::size_t lb_loc = grid_partition.get_local_block(neigh_edge.lb_glob_end);
                            auto bounds = this->get_recv_index_bounds(neigh_edge.edge_vec);
                            std::size_t copy_size = bounds.size(0)*cell_elem_size;
                            for (int k = bounds.min(2); k < bounds.max(2); ++k)
                            {
                                for (int j = bounds.min(1); j < bounds.max(1); ++j)
                                {
                                    std::copy(
                                        (char*)(&(recv_buf_loc[recv_offset_loc])),
                                        (char*)(&(recv_buf_loc[recv_offset_loc]))+copy_size,
                                        (char*)(&array.unwrap_idx(0,bounds.min(0),j,k,lb_loc,maj)));
                                        
                                    recv_offset_loc += copy_size;
                                }
                            }
                        }
                    }
                }
            }
            
            const partition::block_partition_t& get_partition(void) const { return grid_partition; }
        private:
            partition::block_partition_t grid_partition;
            bool is3d;
            coord_t coord_system;
            ctrs::array<dtype,  3> dx;
            ctrs::array<int, 3> num_blocks;
            ctrs::array<int, 3> cells_in_block;
            ctrs::array<int, 3> exchange_cells;
            bound_box_t<dtype,  3> bounds;
            std::vector<bound_box_t<dtype, 3>> block_boxes;
            size_t total_blocks;
            par_group_t* grid_group;
            std::vector<neighbor_relationship_t> neighbors;
            std::vector<ctrs::array<bool, 6>> block_is_domain_boundary;
            std::vector<std::size_t> send_size_elems;
            std::vector<std::size_t> recv_size_elems;
            std::vector<std::vector<char>> send_bufs;
            std::vector<std::vector<char>> recv_bufs;
            std::vector<request_t> requests;
            std::vector<status_t>  statuses;
    };
    
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
        
        template <typename numeric_t> grid_array& operator *= (const numeric_t& rhs)
        requires std::floating_point<numeric_t>
        {
            for (std::size_t i = 0; i < data.size(); ++i) data[i] *= rhs;
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