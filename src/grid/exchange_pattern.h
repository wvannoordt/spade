#pragma once

#include <vector>
#include <tuple>

#include "core/parallel.h"
#include "core/ctrs.h"
#include "core/md_loop.h"
#include "core/poly_vec.h"

#include "grid/grid.h"
#include "grid/partition.h"
#include "grid/transactions.h"
#include "grid/get_transaction.h"

namespace spade::grid
{
    template <typename grid_t>
    struct grid_exchange_config_t
    {
        utils::poly_vec_t<grid_rect_copy_t, patch_fill_t<grid_t::dim()>> send_data;
        utils::poly_vec_t<grid_rect_copy_t, patch_fill_t<grid_t::dim()>> recv_data;
        std::vector<std::size_t>            num_send_elems; //counts the number of CELLS sent or recieved
        std::vector<std::size_t>            num_recv_elems;
        
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
            const auto view = arr.view();
            send_data.foreach([&](const auto& transaction)
            {
                auto& buf    = buffers[transaction.receiver()];
                auto& offset = offsets[transaction.receiver()];
                offset += transaction.insert(view, &buf[offset]);
            });
        }
        
        template <typename array_t>
        void unpack_from(array_t& arr, const std::vector<std::vector<char>>& buffers) const
        {
            std::vector<std::size_t> offsets(num_recv_elems.size(), 0);
            auto view = arr.view();
            recv_data.foreach([&](const auto& transaction)
            {
                const auto& buf = buffers[transaction.sender()];
                auto& offset    = offsets[transaction.sender()];
                offset += transaction.extract(view, &buf[offset]);
            });
        }
        
        template <typename p2p_transaction_t>
        void add_send(const p2p_transaction_t& elem)
        {
            if (elem.sender() == my_rank)
            {
                send_data.add(elem);
                num_send_elems[elem.receiver()] += elem.send_volume();
            }
        }
        
        template <typename p2p_transaction_t>
        void add_recv(const p2p_transaction_t& elem)
        {
            if (elem.receiver() == my_rank)
            {
                recv_data.add(elem);
                num_recv_elems[elem.sender()] += elem.recv_volume();
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
        
        for (auto lb = utils::tag[partition::global](0); lb.value < grid.get_num_global_blocks(); ++lb.value)
        {
            const auto& neighs = grid.get_blocks().get_neighs(lb.value);
            for (const auto& e:neighs)
            {
                bool ignore_from_periodic = false;
                const auto& is_dom_bdy = grid.is_domain_boundary(lb);
                for (int d = 0; d < grid_t::dim(); ++d)
                {
                    bool loc = !is_periodic[d] && ((is_dom_bdy.min(d) && e.edge[d] == -1) || (is_dom_bdy.max(d) && e.edge[d] == 1));
                    ignore_from_periodic = ignore_from_periodic || loc;
                }
                if (!ignore_from_periodic)
                {
                    const auto transaction = get_transaction(grid, lb.value, e);
                    
                    // Note that we add as both and the internal logic inside
                    // these calls will handle the rank checking
                    // if (transaction.reducible())
                    // {
                        // output0.add_send(transaction.reduce());
                        // output0.add_recv(transaction.reduce());
                    // }
                    // else
                    // {
                        output0.add_send(transaction);
                        output0.add_recv(transaction);
                    // }
                }
            }
        }
        
        return array_exchange_t{output0, output1};
    }
}