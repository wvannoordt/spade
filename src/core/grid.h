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

namespace cvdf::grid
{    
    template <class T> concept multiblock_grid = requires(T t, size_t i, size_t j, size_t k, size_t lb)
    {
        // todo: write this
        { t.node_coords(i, j, k, lb) } -> ctrs::vec_nd<3, typename T::dtype>;
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
        t.get_grid();
        t.unwrap_idx(a, i, j, k, lb, b);
    };
    
    enum exchange_inclusion_e
    {
        exclude_exchanges=0,
        include_exchanges=1
    };
    
    template <multiblock_grid grid_t> auto create_grid_dims(const grid_t& grid, const array_center_e& center_type)
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
    
    template <coords::coordinate_system coord_t, parallel::parallel_group par_group_t> class cartesian_grid_t
    {
        public:
            typedef coord_t::coord_type dtype;
            cartesian_grid_t(
                const ctrs::array<std::size_t, cvdf_dim>& num_blocks_in,
                const ctrs::array<std::size_t, cvdf_dim>& cells_in_block_in,
                const ctrs::array<std::size_t, cvdf_dim>& exchange_cells_in,
                const bound_box_t<dtype,  cvdf_dim>& bounds_in,
                const coord_t& coord_system_in,
                par_group_t& group_in)
            {
                //Initialize class members
                this->grid_group = &group_in;
                coord_system = coord_system_in;
                dx = 1.0;
                num_blocks = 1;
                cells_in_block = 1;
                total_blocks = 1;
                is3d = (cvdf_dim==3);
                bounds.min(2) = 0.0;
                bounds.max(2) = 1.0;
                
                //copy constructor arguments, compute domain bounding blox
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
                
                //partition the domain
                grid_partition = partition::block_partition_t(num_blocks, &group_in);
                
                //compute individual block bounding boxes
                block_boxes.resize(grid_partition.get_num_local_blocks());
                for (auto lb: range(0,grid_partition.get_num_local_blocks()))
                {
                    auto& box = block_boxes[lb[0]];
                    ctrs::array<std::size_t, 3> glob_block_idx = ctrs::expand_index(grid_partition.get_global_block(lb[0]), num_blocks);
                    static_for<0,3>([&](auto i)
                    {
                        box.min(i.value) = bounds.min(i.value) + (glob_block_idx[i.value]+0)*bounds.size(i.value)/num_blocks[i.value];
                        box.max(i.value) = bounds.min(i.value) + (glob_block_idx[i.value]+1)*bounds.size(i.value)/num_blocks[i.value];
                    });
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
                    std::size_t lb_glob = lb[0];
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
            
            _finline_ ctrs::array<dtype, 3> node_coords(const int& i, const int& j, const int& k, const int& lb) const
            {
                ctrs::array<dtype, 3> output(0, 0, 0);
                ctrs::array<int, 3> ijk(i, j, k);
                auto& box = block_boxes[lb];
                for (size_t idir = 0; idir < cvdf_dim; idir++)
                {
                    output[idir] = box.min(idir) + ijk[idir]*dx[idir];
                }
                return coord_system.map(output);
            }
            
            _finline_ ctrs::array<dtype, 3> cell_coords(const int& i, const int& j, const int& k, const int& lb) const
            {
                ctrs::array<dtype, 3> output(0, 0, 0);
                ctrs::array<int, 3> ijk(i, j, k);
                auto& box = block_boxes[lb];
                for (size_t idir = 0; idir < cvdf_dim; idir++)
                {
                    output[idir] = box.min(idir) + (ijk[idir]+0.5)*dx[idir];
                }
                return coord_system.map(output);
            }
            
            const par_group_t& group(void) const
            {
                return *grid_group;
            }
            
            const coord_t& coord_sys(void) const {return coord_system;}
            
            md_range_t<int,4> get_range(const array_center_e& centering_in, const exchange_inclusion_e& do_guards=exclude_exchanges) const
            {
                int iexchg = 0;
                if (do_guards==include_exchanges) iexchg = 1;
                switch (centering_in)
                {
                    case cell_centered:
                    {
                        return md_range_t<int,4>(
                            -iexchg*exchange_cells[0],cells_in_block[0]+iexchg*exchange_cells[0],
                            -iexchg*exchange_cells[1],cells_in_block[1]+iexchg*exchange_cells[1],
                            -iexchg*exchange_cells[2],cells_in_block[2]+iexchg*exchange_cells[2],
                            0,grid_partition.get_num_local_blocks());
                    }
                    case node_centered:
                    {
                        return md_range_t<int,4>(
                            -iexchg*exchange_cells[0],1+cells_in_block[0]+iexchg*exchange_cells[0],
                            -iexchg*exchange_cells[1],1+cells_in_block[1]+iexchg*exchange_cells[1],
                            -iexchg*exchange_cells[2],1+cells_in_block[2]+iexchg*exchange_cells[2],
                            0,grid_partition.get_num_local_blocks());
                    }
                    default: return md_range_t<int,4>(0,0,0,0,0,0,0,0);
                }
            }
            
            std::size_t get_num_blocks(const std::size_t& i)   const {return num_blocks[i];}
            std::size_t  get_num_local_blocks(void) const {return grid_partition.get_num_local_blocks();}
            std::size_t get_num_global_blocks(void) const {return num_blocks[0]*num_blocks[1]*num_blocks[2];}
            std::size_t get_num_cells(const std::size_t& i)    const {return cells_in_block[i];}
            std::size_t get_num_exchange(const std::size_t& i) const {return exchange_cells[i];}
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
                for (int i = 0; i < cvdf_dim; i++)
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
                for (int i = 0; i < cvdf_dim; i++)
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
                const std::size_t num_neighs = static_math::pow<3,cvdf_dim>::value;
                
                typedef typename array_t::value_type array_data_t;
                
                std::vector<std::size_t> num_send_bytes(grid_group->size(), 0);
                std::vector<std::size_t> num_recv_bytes(grid_group->size(), 0);
                
                std::size_t cell_elem_size = sizeof(array_data_t)*array.get_minor_dims().total_size();
                
                for (auto p:range(0, grid_group->size()))
                {
                    const std::size_t total_send_buf_size = cell_elem_size*send_size_elems[p[0]];
                    const std::size_t total_recv_buf_size = cell_elem_size*recv_size_elems[p[0]];
                    if (send_bufs[p[0]].size() < total_send_buf_size) send_bufs[p[0]].resize(total_send_buf_size);
                    if (recv_bufs[p[0]].size() < total_recv_buf_size) recv_bufs[p[0]].resize(total_recv_buf_size);
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
                                        (char*)(&array.unwrap_idx(0,bounds.min(0),j,k,lb_loc,maj[0])),
                                        (char*)(&array.unwrap_idx(0,bounds.max(0),j,k,lb_loc,maj[0])),
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
                                        (char*)(&array.unwrap_idx(0,bounds.min(0),j,k,lb_loc,maj[0])));
                                        
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
            ctrs::array<std::size_t, 3> num_blocks;
            ctrs::array<std::size_t, 3> cells_in_block;
            ctrs::array<std::size_t, 3> exchange_cells;
            bound_box_t<dtype,  3> bounds;
            std::vector<bound_box_t<dtype, 3>> block_boxes;
            size_t total_blocks;
            par_group_t* grid_group;
            std::vector<neighbor_relationship_t> neighbors;
            
            std::vector<std::size_t> send_size_elems;
            std::vector<std::size_t> recv_size_elems;
            std::vector<std::vector<char>> send_bufs;
            std::vector<std::vector<char>> recv_bufs;
            std::vector<request_t> requests;
            std::vector<status_t>  statuses;
    };
    
    template <
        multiblock_grid grid_t,
        typename ar_data_t,
        dims::grid_array_dimension minor_dim_t,
        dims::grid_array_dimension major_dim_t,
        const array_center_e centering=cell_centered,
        array_container::grid_data_container container_t=std::vector<ar_data_t>>
    struct grid_array
    {
        grid_array(
            const grid_t& grid_in,
            const ar_data_t& fill_elem,
            const minor_dim_t& minor_dims_in,
            const major_dim_t& major_dims_in)
        {
            grid = &grid_in;
            minor_dims = minor_dims_in;
            major_dims = major_dims_in;
            grid_dims = create_grid_dims(grid_in, this->centering_type());
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
        
        typedef grid_t grid_type;
        typedef ar_data_t value_type;
        
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
        
        template <typename idx_t> _finline_ std::size_t offst_r(const std::size_t& lev, const idx_t& idx) const
        {
            static_assert(std::is_integral<idx_t>::value, "integral type required for array indexing");
            return idx_coeffs[lev]*idx;
        }
        
        template <typename idx_t, typename... idxs_t> _finline_ std::size_t offst_r(const std::size_t& lev, const idx_t& idx, idxs_t... idxs) const
        {
            static_assert(std::is_integral<idx_t>::value, "integral type required for array indexing");
            return idx_coeffs[lev]*idx + offst_r(lev+1, idxs...);
        }
        
        template <typename... idxs_t>
        _finline_ ar_data_t& operator() (idxs_t... idxs)
        {
            static_assert(sizeof...(idxs_t)==total_idx_rank(), "wrong number of indices passed to indexing operator");
            return data[offset + offst_r(0, idxs...)];
        }
        
        template <typename... idxs_t>
        _finline_ const ar_data_t& operator() (idxs_t... idxs) const
        {
            static_assert(sizeof...(idxs_t)==total_idx_rank(), "wrong number of indices passed to indexing operator");
            return data[offset + offst_r(0, idxs...)];
        }
        
        template <typename i0_t, typename i1_t, typename i2_t, typename i3_t, typename i4_t, typename i5_t>
        _finline_ ar_data_t& unwrap_idx(const i0_t& i0, const i1_t& i1, const i2_t& i2, const i3_t& i3, const i4_t& i4, const i5_t& i5)
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
        
        template <typename i0_t, typename i1_t, typename i2_t, typename i3_t, typename i4_t, typename i5_t>
        _finline_ const ar_data_t& unwrap_idx(const i0_t& i0, const i1_t& i1, const i2_t& i2, const i3_t& i3, const i4_t& i4, const i5_t& i5) const
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