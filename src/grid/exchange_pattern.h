#pragma once

#include <vector>
#include <tuple>

#include "grid/grid.h"
#include "core/parallel.h"
#include "core/ctrs.h"
#include "grid/partition.h"
#include "core/md_loop.h"

namespace spade::grid
{
    struct grid_rect_copy_t
    {
        bound_box_t<int, 4> source, dest;
        int rank_recv, rank_send;
        std::size_t volume() const { return source.volume(); }
    };

    template <typename grid_t>
    struct grid_exchange_config_t
    {
        std::vector<grid_rect_copy_t> send_data;
        std::vector<grid_rect_copy_t> recv_data;
        std::vector<std::size_t>      num_send_elems; //counts the number of CELLS sent or recieved
        std::vector<std::size_t>      num_recv_elems;
        
        int my_rank;
        const grid_t& grid;
        grid_exchange_config_t(const grid_t& grid_in) : grid{grid_in}, num_send_elems{0}, num_recv_elems{0}
        {
            auto np = grid.group().size();
            num_send_elems.resize(np);
            num_recv_elems.resize(np);
            my_rank = grid.group().rank();
        }
        
        template <typename array_t>
        std::size_t get_num_recv_bytes(const array_t& arr) const
        {
            const std::size_t fundamentals_per_cell = arr.cell_element_size();
            const std::size_t bytes_per_fundamental = sizeof(typename array_t::fundamental_type);
            return num_recv_elems[arr.get_grid().group().rank()] * fundamentals_per_cell * bytes_per_fundamental;
        }
        
        template <typename array_t>
        std::size_t get_num_send_bytes(const array_t& arr) const
        {
            const std::size_t fundamentals_per_cell = arr.cell_element_size();
            const std::size_t bytes_per_fundamental = sizeof(typename array_t::fundamental_type);
            return num_send_elems[arr.get_grid().group().rank()] * fundamentals_per_cell * bytes_per_fundamental;
        }
        
        template <typename array_t>
        void pack_to(const array_t& arr, std::vector<std::vector<char>>& buffers) const
        {
            std::vector<std::size_t> offsets(num_send_elems.size(), 0);
            using element_type = typename array_t::alias_type;
            element_type elem;
            char* raw_e_addr = (char*)&elem;
            std::size_t elem_size = sizeof(element_type);
            int count = 0;
            int count2 = 0;
            for (const auto& transaction: send_data)
            {
                auto& buf    = buffers[transaction.rank_recv];
                auto& offset = offsets[transaction.rank_recv];
                grid::cell_idx_t idx;
                algs::md_loop(idx, transaction.source, [&](const auto& icell)
                {
                    elem = arr.get_elem(idx);
                    std::copy(raw_e_addr, raw_e_addr + elem_size, &buf[offset]);
                    offset += elem_size;
                });
            }
        }
        
        template <typename array_t>
        void unpack_from(array_t& arr, const std::vector<std::vector<char>>& buffers) const
        {
            std::vector<std::size_t> offsets(num_recv_elems.size(), 0);
            using element_type = typename array_t::alias_type;
            element_type elem;
            char* raw_e_addr = (char*)&elem;
            std::size_t elem_size = sizeof(element_type);
            for (const auto& transaction: recv_data)
            {
                const auto& buf = buffers[transaction.rank_send];
                auto& offset    = offsets[transaction.rank_send];
                grid::cell_idx_t idx;
                algs::md_loop(idx, transaction.dest, [&](const auto& icell)
                {
                    std::copy(&buf[offset], &buf[offset]+elem_size, raw_e_addr);
                    arr.set_elem(idx, elem);
                    offset += elem_size;
                });
            }
        }
        
        void add_send(const grid_rect_copy_t& elem)
        {
            if (elem.rank_send == my_rank)
            {
                send_data.push_back(elem);
                num_send_elems[elem.rank_recv] += elem.volume();
            }
        }
        
        void add_recv(const grid_rect_copy_t& elem)
        {
            if (elem.rank_recv == my_rank)
            {
                recv_data.push_back(elem);
                num_recv_elems[elem.rank_send] += elem.volume();
            }
        }
    };
    
    
    template <typename par_group_t>
    struct exchange_handle_t
    {
        par_group_t& group;
        std::vector<std::vector<char>> send_bufs;
        std::vector<std::vector<char>> recv_bufs;
        std::vector<request_t> requests;
        std::vector<status_t>  statuses;
        
        exchange_handle_t(par_group_t& group_in) : group{group_in}
        {
            send_bufs.resize(group.size());
            recv_bufs.resize(group.size());
            
            requests.resize(group.size());
            statuses.resize(group.size());
        }
        
        template <typename array_t, typename grid_t>
        requires(std::same_as<typename array_t::grid_type, grid_t>)
        void assure_buffer_size(const array_t& arr, const grid_exchange_config_t<grid_t>& config)
        {
            const std::size_t elem_size = sizeof(typename array_t::alias_type);
            for (int p = 0; p < group.size(); ++p)
            {
                std::size_t num_send_bytes = config.num_send_elems[p]*elem_size;
                std::size_t num_recv_bytes = config.num_recv_elems[p]*elem_size;
                
                if (send_bufs[p].size() < num_send_bytes) send_bufs[p].resize(num_send_bytes);
                if (recv_bufs[p].size() < num_recv_bytes) recv_bufs[p].resize(num_recv_bytes);
            }
        }
    };
    
    template <typename array_t, typename config_t, typename handle_t>
    struct handle_sentinel_t
    {
        array_t& array;
        config_t& config;
        handle_t& handle;
        
        void await()
        {
            handle.group.await_all(handle.statuses.size(), handle.requests.data(), handle.statuses.data());
            config.unpack_from(array, handle.recv_bufs);
        }
    };
    
    template <typename grid_t, typename par_group_t>
    struct array_exchange_t
    {
        grid_exchange_config_t<grid_t> config;
        exchange_handle_t<par_group_t> handle;
        
        template <typename arr_t>
        [[nodiscard("exchange sentinel must be used to complete exchange")]]
        auto async_exchange(arr_t& array)
        {
            handle.assure_buffer_size(array, config);
            config.pack_to(array, handle.send_bufs);
            for (std::size_t p = 0; p < handle.group.size(); ++p)
            {
                request_t r = handle.group.async_recv(&handle.recv_bufs[p][0], handle.recv_bufs[p].size(), p);
                handle.requests[p] = r;
            }
            for (std::size_t p = 0; p < handle.group.size(); ++p)
            {
                handle.group.sync_send(&handle.send_bufs[p][0], handle.send_bufs[p].size(), p);
            }
            return handle_sentinel_t{array, config, handle};
        }
        
        template <typename arr_t>
        void exchange(arr_t& array)
        {
            auto handle = this->async_exchange(array);
            handle.await();
        }
    };
    
    
    template <typename grid_t, typename group_t>
    auto create_exchange(const grid_t& grid, group_t& group, const ctrs::array<bool, grid_t::dim()>& is_periodic)
    {
        grid_exchange_config_t output0(grid);
        exchange_handle_t      output1(group);
        
        const auto get_transaction = [&](const auto& relation)
        {
            grid_rect_copy_t output;
            
            auto lb_ini  = utils::tag[partition::global](relation.lb_ini);
            auto lb_term = utils::tag[partition::global](relation.lb_term);
            
            output.rank_send = grid.get_partition().get_rank(lb_ini);
            output.rank_recv = grid.get_partition().get_rank(lb_term);
            
            output.source.min(3) = grid.get_partition().to_local(lb_ini).value;
            output.source.max(3) = output.source.min(3)+1;

            output.dest.min(3)   = grid.get_partition().to_local(lb_term).value;
            output.dest.max(3)   = output.dest.min(3)+1;
            
            for (int d = 0; d < 3; ++d)
            {
                const int sign = relation.edge[d];
                const int nx   = grid.get_num_cells(d);
                const int ng   = grid.get_num_exchange(d);
                switch (sign)
                {
                    case -1:
                    {
                        output.source.min(d) = 0;
                        output.source.max(d) = ng;
                        output.dest.min(d)   = nx;
                        output.dest.max(d)   = nx+ng;
                        break;
                    }
                    case  0:
                    {
                        output.source.min(d) = 0;
                        output.source.max(d) = nx;
                        output.dest.min(d)   = 0;
                        output.dest.max(d)   = nx;
                        break;
                    }
                    case  1:
                    {
                        output.source.min(d) = nx-ng;
                        output.source.max(d) = nx;
                        output.dest.min(d)   = -ng;
                        output.dest.max(d)   = 0;
                        break;
                    }
                }
            }
            return output;
        };
        
        for (auto lb = utils::tag[partition::global](0); lb.value < grid.get_num_global_blocks(); ++lb.value)
        {
            const auto neighs = grid.get_blocks().get_neighs(lb.value);
            for (const auto& e:neighs)
            {
                grid_rect_copy_t transaction = get_transaction(e);
                
                // Note that we add as both and the internal logic inside
                // these calls will handle the rank checking
                output0.add_send(transaction);
                output0.add_recv(transaction);
            }
        }
        
        return array_exchange_t{output0, output1};
    }
}