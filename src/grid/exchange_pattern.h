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
        const grid_t* grid;
        grid_exchange_config_t(const grid_t& grid_in) : grid{&grid_in}, num_send_elems{0}, num_recv_elems{0}
        {
            auto np = grid->group().size();
            num_send_elems.resize(np);
            num_recv_elems.resize(np);
            my_rank = grid->group().rank();
        }
        grid_exchange_config_t(){} 
        
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
            const auto img = arr.image();
            send_data.foreach([&](const auto& transaction)
            {
                auto& buf    = buffers[transaction.receiver()];
                auto& offset = offsets[transaction.receiver()];
                offset += transaction.insert(img, &buf[offset]);
            });
        }
        
        template <typename array_t>
        void unpack_from(array_t& arr, const std::vector<std::vector<char>>& buffers) const
        {
            std::vector<std::size_t> offsets(num_recv_elems.size(), 0);
            auto img = arr.image();
            recv_data.foreach([&](const auto& transaction)
            {
                const auto& buf = buffers[transaction.sender()];
                auto& offset    = offsets[transaction.sender()];
                offset += transaction.extract(img, &buf[offset]);
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
        const par_group_t* group;
        std::vector<std::vector<char>> send_bufs;
        std::vector<std::vector<char>> recv_bufs;
        std::vector<request_t> requests;
        std::vector<status_t>  statuses;
        
        exchange_handle_t(const par_group_t& group_in) : group{&group_in}
        {
            send_bufs.resize(group->size());
            recv_bufs.resize(group->size());
            
            requests.resize(group->size());
            statuses.resize(group->size());
        }
        
        exchange_handle_t(){}
        
        template <typename array_t, typename grid_t>
        requires(std::same_as<typename array_t::grid_type, grid_t>)
        void assure_buffer_size(const array_t& arr, const grid_exchange_config_t<grid_t>& config)
        {
            const std::size_t elem_size = sizeof(typename array_t::alias_type);
            for (int p = 0; p < group->size(); ++p)
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
            handle.group->await_all(handle.statuses.size(), handle.requests.data(), handle.statuses.data());
            config.unpack_from(array, handle.recv_bufs);
        }
    };
    
    template <typename grid_t>
    struct gpu_exch_t
    {
        constexpr static int interp_size = static_math::pow<2, grid_t::dim()>::value;
        device::shared_vector<std::size_t> inj_sends;
        device::shared_vector<std::size_t> inj_recvs;
        device::shared_vector<std::size_t> int_sends;
        device::shared_vector<std::size_t> int_recvs;
        bool configd = false;
    };
    
    template <typename grid_t, typename par_group_t>
    struct array_exchange_t// : public grid_t::dependent_type
    {
        grid_exchange_config_t<grid_t> config;
        exchange_handle_t<par_group_t> handle;
        gpu_exch_t<grid_t>             device_exch;
        
        array_exchange_t(const grid_exchange_config_t<grid_t>& cfgin, const exchange_handle_t<par_group_t>& hdlin)
        : config{cfgin}, handle{hdlin}
        {

        }
        
        array_exchange_t(){}
        
        template <typename arr_t>
        void config_gpu_exch(const arr_t& arr)
        {
            typename arr_t::var_idx_t i0 = 0;
            
            const auto aimg = arr.image();
            
            const auto get_base_offst = [&](const grid::cell_idx_t& icell)
            {
                return aimg.map.compute_offset(i0, icell);
            };
            
            const auto& injecs = config.send_data.data;
            const auto& interp = config.send_data.sub.data;
            for (const auto& injc: injecs)
            {
                const auto& snds = injc.source;
                const auto& recv = injc.dest;
                
                const auto l0 = [&](const auto& arr_i)
                {
                    const auto loc_oft = get_base_offst(grid::cell_idx_t(arr_i[0], arr_i[1], arr_i[2], arr_i[3]));
                    device_exch.inj_sends.push_back(loc_oft);
                };
                const auto l1 = [&](const auto& arr_i)
                {
                    const auto loc_oft = get_base_offst(grid::cell_idx_t(arr_i[0], arr_i[1], arr_i[2], arr_i[3]));
                    device_exch.inj_recvs.push_back(loc_oft);
                };
                dispatch::execute(snds, l0, device::cpu);
                dispatch::execute(recv, l1, device::cpu);
            }
            for (const auto& intp: interp)
            {
                const auto& recvs = intp.patches.dest;
                const auto l0 = [&](const auto& arr_i)
                {
                    grid::cell_idx_t gidx(arr_i[0], arr_i[1], arr_i[2], arr_i[3]);
                    int di = gidx.i() - intp.patches.dest.min(0);
                    int dj = gidx.j() - intp.patches.dest.min(1);
                    int dk = gidx.k() - intp.patches.dest.min(2);
                    algs::static_for<0,intp.max_size>([&](const auto iii)
                    {
                        constexpr int iv = iii.value;
                        constexpr int d0 = (iv >> 0)&1;
                        constexpr int d1 = (iv >> 1)&1;
                        constexpr int d2 = (iv >> 2)&1;
                        
                        int i0 = di<<(intp.i_coeff[0]+1);
                        int j0 = dj<<(intp.i_coeff[1]+1);
                        int k0 = dk<<(intp.i_coeff[2]+1);
                        i0 = i0 >> 1;
                        j0 = j0 >> 1;
                        k0 = k0 >> 1;
                        
                        
                        const int donor_di = i0 + d0*intp.i_incr[0];
                        const int donor_dj = j0 + d1*intp.i_incr[1];
                        const int donor_dk = k0 + d2*intp.i_incr[2];
                        
                        grid::cell_idx_t donor_idx;
                        donor_idx.lb() = intp.patches.source.min(3);
                        donor_idx.i()  = intp.patches.source.min(0) + donor_di;
                        donor_idx.j()  = intp.patches.source.min(1) + donor_dj;
                        donor_idx.k()  = intp.patches.source.min(2) + donor_dk;
                        device_exch.int_sends.push_back(get_base_offst(donor_idx));
                    });
                    const auto loc_oft = get_base_offst(gidx);
                    device_exch.int_recvs.push_back(loc_oft);
                };
                dispatch::execute(recvs, l0, device::cpu);
            }
            
            device_exch.inj_sends.transfer();
            device_exch.inj_recvs.transfer();
            device_exch.int_sends.transfer();
            device_exch.int_recvs.transfer();
            device_exch.configd = true;
        }
        
        template <typename arr_t>
        [[nodiscard("exchange sentinel must be used to complete exchange")]]
        auto async_exchange(arr_t& array)
        {
            handle.assure_buffer_size(array, config);
            config.pack_to(array, handle.send_bufs);
            for (std::size_t p = 0; p < handle.group->size(); ++p)
            {
                request_t r = handle.group->async_recv(&handle.recv_bufs[p][0], handle.recv_bufs[p].size(), p);
                handle.requests[p] = r;
            }
            for (std::size_t p = 0; p < handle.group->size(); ++p)
            {
                handle.group->sync_send(&handle.send_bufs[p][0], handle.send_bufs[p].size(), p);
            }
            return handle_sentinel_t{array, config, handle};
        }
        
        template <typename arr_t>
        requires device::is_cpu<typename arr_t::device_type>
        void exchange(arr_t& array)
        {
            auto handle = this->async_exchange(array);
            handle.await();
        }
        
        template <typename arr_t>
        requires device::is_gpu<typename arr_t::device_type>
        void exchange(const arr_t& src_arr, arr_t& dst_arr)
        {
            if (!device_exch.configd)
            {
                config_gpu_exch(src_arr);
                device_exch.configd = true;
            }
            using data_t = typename arr_t::fundamental_type;
            using alias_t = typename arr_t::alias_type;
            const auto src_img = utils::make_vec_image(src_arr.data);
            auto       dst_img = utils::make_vec_image(dst_arr.data);
            cell_idx_t root_idx(0,0,0,0);
            
            auto img = dst_arr.image();
            
            constexpr std::size_t nv = [&]()
            {
                if constexpr (ctrs::basic_array<alias_t>) return alias_t::size();
                else                                      return 1;
            }();
            
            const std::size_t voffst = [&]()
            {
                if constexpr (ctrs::basic_array<alias_t>) return img.map.compute_offset(1, root_idx) - img.map.compute_offset(0, root_idx);
                else                                      return 0;
            }();
            
            const auto i_inj_s = utils::make_vec_image(device_exch.inj_sends.data(src_arr.device()));
            const auto i_inj_r = utils::make_vec_image(device_exch.inj_recvs.data(src_arr.device()));
            const auto i_int_s = utils::make_vec_image(device_exch.int_sends.data(src_arr.device()));
            const auto i_int_r = utils::make_vec_image(device_exch.int_recvs.data(src_arr.device()));
            
            dispatch::ranges::linear_range_t injections    (std::size_t(0), device_exch.inj_recvs.data(src_arr.device()).size(), src_arr.device());
            dispatch::ranges::linear_range_t interpolations(std::size_t(0), device_exch.int_recvs.data(src_arr.device()).size(), src_arr.device());
            
            auto inj_load = _sp_lambda (const std::size_t& idx) mutable
            {
                for (std::size_t v = 0; v < nv; ++v)
                {
                    dst_img[i_inj_r[idx] + v*voffst] = src_img[i_inj_s[idx] + v*voffst];
                }
            };
            
            const int ndonor = device_exch.interp_size;
            const data_t coeff = data_t(1.0)/ndonor;
            auto int_load = _sp_lambda (const std::size_t& idx) mutable
            {
                for (std::size_t v = 0; v < nv; ++v)
                {
                    dst_img[i_int_r[idx] + v*voffst] = data_t(0);
                    for (std::size_t id = 0; id < ndonor; ++id)
                    {
                        dst_img[i_int_r[idx] + v*voffst] += coeff*src_img[i_int_s[id + ndonor*idx] + v*voffst];
                    }
                }
            };
            
            dispatch::execute(interpolations, int_load);
            dispatch::execute(injections,     inj_load);
        }
        
        template <typename arr_t>
        requires device::is_gpu<typename arr_t::device_type>
        void exchange(arr_t& arr)
        {
            exchange(arr, arr);
        }
    };
    
    
    template <typename grid_t, typename group_t>
    auto create_exchange(const grid_t& grid, const group_t& group, const ctrs::array<bool, grid_t::dim()>& is_periodic)
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
                    if (transaction.reducible())
                    {
                        output0.add_send(transaction.reduce());
                        output0.add_recv(transaction.reduce());
                    }
                    else
                    {
                        output0.add_send(transaction);
                        output0.add_recv(transaction);
                    }
                }
            }
        }
        
        return array_exchange_t(output0, output1);
    }
}