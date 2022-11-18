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
    const static static_math::int_const_t<3> dim3;
    const static static_math::int_const_t<2> dim2;
    
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
    
    template <class T, const array_centering ct> concept has_centering_type = (T::centering_type() == ct);
    
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
    
    template <multiblock_grid grid_t, const array_centering centering>
    requires (centering == cell_centered)
    static auto create_grid_dims(const grid_t& grid)
    {
        return dims::dynamic_dims<4>(
            grid.get_num_cells(0)+grid.get_num_exchange(0)*2,
            grid.get_num_cells(1)+grid.get_num_exchange(1)*2,
            grid.get_num_cells(2)+grid.get_num_exchange(2)*2,
            grid.get_num_local_blocks());
    }
    
    template <multiblock_grid grid_t, const array_centering centering>
    requires (centering == node_centered)
    static auto create_grid_dims(const grid_t& grid)
    {
        return dims::dynamic_dims<4>(
            1+grid.get_num_cells(0)+grid.get_num_exchange(0)*2,
            1+grid.get_num_cells(1)+grid.get_num_exchange(1)*2,
            1+grid.get_num_cells(2)+grid.get_num_exchange(2)*2,
            grid.get_num_local_blocks());
    }
    
    template <multiblock_grid grid_t, const array_centering centering>
    requires (centering == face_centered)
    static auto create_grid_dims(const grid_t& grid)
    {
        return dims::dynamic_dims<5>(
            1+grid.get_num_cells(0)+grid.get_num_exchange(0)*2,
            1+grid.get_num_cells(1)+grid.get_num_exchange(1)*2,
            1+grid.get_num_cells(2)+grid.get_num_exchange(2)*2,
            grid.dim(),
            grid.get_num_local_blocks());
    }
    
    struct neighbor_relationship_t
    {
        ctrs::array<int, 3> edge_vec;
        int rank_end, rank_start;
        int lb_glob_end, lb_glob_start;
    };
    
    template
    <
        coords::coordinate_system coord_t,
        parallel::parallel_group par_group_t,
        const std::size_t grid_dim
    >
    class cartesian_grid_t
    {
        public:
            typedef coord_t::coord_type dtype;
            typedef coord_t::coord_type coord_type;
            typedef coord_t coord_sys_type;
            typedef ctrs::array<coord_type, 3> coord_point_type;
            cartesian_grid_t(
                const ctrs::array<int,   grid_dim>& num_blocks_in,
                const ctrs::array<int,   grid_dim>& cells_in_block_in,
                const ctrs::array<int,   grid_dim>& exchange_cells_in,
                const bound_box_t<dtype, grid_dim>& bounds_in,
                const coord_t& coord_system_in,
                par_group_t& group_in
                )
            {
                //Initialize class members
                this->grid_group = &group_in;
                coord_system = coord_system_in;
                dx = 1.0;
                num_blocks = 1;
                cells_in_block = 1;
                exchange_cells = 0;
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
                    ctrs::array<int, 3> glob_block_idx = ctrs::expand_index(grid_partition.get_global_block(lb), num_blocks);
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
                    for (auto d: range(0,dim()))
                    {
                        data.min(d) = (lb_idx[d]==0);
                        data.max(d) = (lb_idx[d]==(num_blocks[d]-1));
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
            
            cartesian_grid_t(){}
            
            constexpr static int dim(void) {return grid_dim;}
            
            template <typename idx_t>_finline_ ctrs::array<dtype, 3> get_coords(const idx_t& i) const
            {
                return coord_system.map(this->get_comp_coords(i));
            }
            
            template <typename idx_t>_finline_ ctrs::array<dtype, 3> get_comp_coords(const idx_t& i) const
            {
                const auto  idx_r = get_index_coord(i);
                const int   lb    = i.lb();
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
            
            md_range_t<int,4> get_range(const array_centering& centering_in, const exchange_inclusion_e& do_guards=exclude_exchanges) const
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
            ctrs::array<int, 3> get_num_blocks(void) const {return num_blocks;}
            std::size_t  get_num_local_blocks(void) const {return grid_partition.get_num_local_blocks();}
            std::size_t get_num_global_blocks(void) const {return num_blocks[0]*num_blocks[1]*num_blocks[2];}
            int get_num_cells(const std::size_t& i)    const {return cells_in_block[i];}
            int get_num_exchange(const std::size_t& i) const {return exchange_cells[i];}
            const bound_box_t<bool, dim()>& is_domain_boundary(const std::size_t& lb_glob) const {return block_is_domain_boundary[lb_glob];}
            bound_box_t<dtype,  3> get_bounds(void) const {return bounds;}
            ctrs::array<dtype,  3> get_dx(void) const {return dx;}
            dtype get_dx(const std::size_t& i) const {return dx[i];}
            bound_box_t<dtype, 3> get_block_box(const std::size_t& lb) const {return block_boxes[lb];}
            bool is_3d(void) const {return this->is3d;}
            
            std::size_t grid_size() const
            {
                return grid_partition.get_num_global_blocks()*cells_in_block[0]*cells_in_block[1]*cells_in_block[2];
            }
            
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
            std::vector<bound_box_t<bool, dim()>> block_is_domain_boundary;
            std::vector<std::size_t> send_size_elems;
            std::vector<std::size_t> recv_size_elems;
            std::vector<std::vector<char>> send_bufs;
            std::vector<std::vector<char>> recv_bufs;
            std::vector<request_t> requests;
            std::vector<status_t>  statuses;
    };
}