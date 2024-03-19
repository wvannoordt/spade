#pragma once

#include <vector>
#include <tuple>

#include "core/parallel.h"
#include "core/ctrs.h"
#include "core/md_loop.h"
#include "core/poly_vec.h"

#include "amr/amr.h"

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
        
        const grid_t* grid0;
        const grid_t* grid1;
        grid_exchange_config_t(const grid_t& grid_in) : grid0{&grid_in}, grid1{&grid_in}, num_send_elems{0}, num_recv_elems{0}
        {
            auto np = grid0->group().size();
            num_send_elems.resize(np);
            num_recv_elems.resize(np);
            my_rank = grid0->group().rank();
        }
        grid_exchange_config_t(const grid_t& grid0_in, const grid_t& grid1_in) : grid0{&grid0_in}, grid1{&grid1_in}, num_send_elems{0}, num_recv_elems{0}
        {
            auto np = grid0->group().size();
            num_send_elems.resize(np);
            num_recv_elems.resize(np);
            my_rank = grid0->group().rank();
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
        
        //Note that these may be relative to different grids
        //Convention: all sends are from one grid, all recieves are from the other
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
        void config_gpu_exch(const arr_t& src, const arr_t& dst)
        {
            typename arr_t::var_idx_t i0 = 0;
            
            const auto src_img = src.image();
            const auto dst_img = dst.image();
            const auto get_base_offst = [&](const auto& img, const grid::cell_idx_t& icell)
            {
                return img.map.compute_offset(i0, icell);
            };
            
            const auto& injecs = config.send_data.data;
            const auto& interp = config.send_data.sub.data;
            for (const auto& injc: injecs)
            {
                const auto& snds = injc.source;
                const auto& recv = injc.dest;
                
                const auto l0 = [&](const auto& arr_i)
                {
                    const auto loc_oft = get_base_offst(src_img, grid::cell_idx_t(arr_i[0], arr_i[1], arr_i[2], arr_i[3]));
                    device_exch.inj_sends.push_back(loc_oft);
                };
                const auto l1 = [&](const auto& arr_i)
                {
                    const auto loc_oft = get_base_offst(dst_img, grid::cell_idx_t(arr_i[0], arr_i[1], arr_i[2], arr_i[3]));
                    device_exch.inj_recvs.push_back(loc_oft);
                };
                dispatch::execute(snds, l0, device::cpu);
                dispatch::execute(recv, l1, device::cpu);
            }
            for (const auto& intp: interp)
            {
                const auto& recvs = intp.patches.dest;
                const auto  l0 = [&](const auto& arr_i)
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
                        device_exch.int_sends.push_back(get_base_offst(src_img, donor_idx));
                    });
                    const auto loc_oft = get_base_offst(dst_img, gidx);
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
        void config_gpu_exch(const arr_t& arr)
        {
            config_gpu_exch(arr, arr);
        }
        
        template <typename arr_t>
        [[nodiscard("exchange sentinel must be used to complete exchange")]]
        auto async_exchange(const arr_t& src_array, arr_t& dst_array)
        {
            handle.assure_buffer_size(src_array, config);
            config.pack_to(src_array, handle.send_bufs);
            for (std::size_t p = 0; p < handle.group->size(); ++p)
            {
                request_t r = handle.group->async_recv(&handle.recv_bufs[p][0], handle.recv_bufs[p].size(), p);
                handle.requests[p] = r;
            }
            for (std::size_t p = 0; p < handle.group->size(); ++p)
            {
                handle.group->sync_send(&handle.send_bufs[p][0], handle.send_bufs[p].size(), p);
            }
            return handle_sentinel_t{dst_array, config, handle};
        }
        
        template <typename arr_t>
        requires device::is_cpu<typename arr_t::device_type>
        void exchange(const arr_t& src, arr_t& dst)
        {
            auto handle = this->async_exchange(src, dst);
            handle.await();
        }
        
        template <typename arr_t>
        requires device::is_cpu<typename arr_t::device_type>
        void exchange(arr_t& arr)
        {
            auto handle = this->async_exchange(arr, arr);
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
            const auto i_inj_r = utils::make_vec_image(device_exch.inj_recvs.data(dst_arr.device()));
            const auto i_int_s = utils::make_vec_image(device_exch.int_sends.data(src_arr.device()));
            const auto i_int_r = utils::make_vec_image(device_exch.int_recvs.data(dst_arr.device()));
            
            dispatch::ranges::linear_range_t injections    (std::size_t(0), device_exch.inj_recvs.data(dst_arr.device()).size(), dst_arr.device());
            dispatch::ranges::linear_range_t interpolations(std::size_t(0), device_exch.int_recvs.data(dst_arr.device()).size(), dst_arr.device());
            
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
    
    
    template <typename array_t, typename group_t>
    auto create_exchange(const array_t& array, const group_t& group, const ctrs::array<bool, array_t::grid_type::dim()>& is_periodic)
    {
        const auto& grid = array.get_grid();
        using grid_t = typename array_t::grid_type;
        const auto  num_exchanges = array.get_num_exchange();
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
                    const auto transaction = get_transaction(num_exchanges, grid, grid, lb.value, e);
                    
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
    
    template <typename arr_t, typename group_t>
    auto create_interface(const arr_t& src_arr, const arr_t& dst_arr, const group_t& group)
    {
        using grid_t = typename arr_t::grid_type;
        static_assert(grid_t::dim() == 3, "can only perform interface generation on 3d grids");
        using required_blocks_type = amr::amr_blocks_t<typename grid_t::coord_type, typename grid_t::array_desig_type>;
        static_assert(std::same_as<typename grid_t::blocks_type, required_blocks_type>, "can only perform interface generation grids with AMR block configs");
        
        const auto& src_grid = src_arr.get_grid();
        const auto& dst_grid = dst_arr.get_grid();
        
        grid_exchange_config_t output0(src_grid, dst_grid);
        exchange_handle_t      output1(group);

        // todo: something more sophisticated than this...
        for (int i = 0; i < grid_t::dim(); ++i)
        {
            // Note: could probably do a better job of exceptions here!
            if (src_grid.get_num_cells(i) != dst_grid.get_num_cells(i))
            {
                std::stringstream ss;
                ss << "attempted to create interface between incompatible grids:";
                ss << "number of grid cells does not match\n";
                ss << "Num. cells, donor: [";
                ss << src_grid.get_num_cells(0) << ", ";
                ss << src_grid.get_num_cells(1) << ", ";
                ss << src_grid.get_num_cells(2) << "]\n";
                ss << "Num. cells, recvr: [";
                ss << dst_grid.get_num_cells(0) << ", ";
                ss << dst_grid.get_num_cells(1) << ", ";
                ss << dst_grid.get_num_cells(2) << "]\n";
                
                throw except::sp_exception(ss.str());
            }
            if (src_arr.get_num_exchange(i) != dst_arr.get_num_exchange(i))
            {
                std::stringstream ss;
                ss << "attempted to create interface between incompatible grids:";
                ss << "number of exchange cells does not match\n";
                ss << "Num. exchange, donor: [";
                ss << src_arr.get_num_exchange(0) << ", ";
                ss << src_arr.get_num_exchange(1) << ", ";
                ss << src_arr.get_num_exchange(2) << "]\n";
                ss << "Num. exchange, recvr: [";
                ss << dst_arr.get_num_exchange(0) << ", ";
                ss << dst_arr.get_num_exchange(1) << ", ";
                ss << dst_arr.get_num_exchange(2) << "]\n";
                throw except::sp_exception(ss.str());
            }
        }
        
        using real_t = typename grid_t::coord_type;
        const real_t tol = 1e-5;
        const auto rt_eq = [&](const real_t a, const real_t b) { return utils::abs(a-b)<tol; };
        int match_dir = -1;
        int src_pm = -1;
        int dst_pm = -1;
        const auto& src_bnd = src_grid.get_bounds();
        const auto& dst_bnd = dst_grid.get_bounds();
        
        // Note: need more checks in here to ensure this doesn't completely break
        for (int i = 0; i < grid_t::dim(); ++i)
        {
            if (rt_eq(src_bnd.min(i), dst_bnd.max(i)))
            {
                match_dir = i;
                src_pm    = 0;
                dst_pm    = 1;
            }
            if (rt_eq(src_bnd.max(i), dst_bnd.min(i)))
            {
                match_dir = i;
                src_pm    = 1;
                dst_pm    = 0;
            }
        }
        if (match_dir < 0)
        {
            throw except::sp_exception("Attempted to create interface, but grids do not align");
        }
        
        // match_dir is the direction in which the interface is normal
        // src_pm = 1 if the interface is on the maximal face in "i" direction relative to the source grid, etc.
        
        const auto mk_lam = [&](const int& idir, const int& pm, const grid_t& grid)
        {
            
            const auto output = [&](const auto& lb)
            {
                const auto& idb = grid.is_domain_boundary(lb);
                return idb(idir, pm);
            };
            return output;
        };
        
        const auto src_check = mk_lam(match_dir, src_pm, src_grid);
        const auto dst_check = mk_lam(match_dir, dst_pm, dst_grid);
        const auto src_lbs = src_grid.select_blocks(src_check, partition::global);
        const auto dst_lbs = dst_grid.select_blocks(dst_check, partition::global);
        
        const auto extended_bbx = [&](const auto& arr, const auto& lb)
        {
            const auto& grd = arr.get_grid();
            auto bbx      = grd.get_bounding_box(lb);
            const auto dx = grd.get_dx(lb);
            for (int d = 0; d < grd.dim(); ++d)
            {
                bbx.min(d) -= arr.get_num_exchange(d)*dx[d];
                bbx.max(d) += arr.get_num_exchange(d)*dx[d];
            }
            return bbx;
        };
        
        for (const auto& src_lb: src_lbs)
        {
            const auto src_bbx = src_grid.get_bounding_box(src_lb);
            for (const auto& dst_lb: dst_lbs)
            {
                const auto dst_bbx_ext = extended_bbx(dst_arr, dst_lb);
                if (src_bbx.intersects(dst_bbx_ext))
                {
                    // src_lb donates to dst_lb
                    // We create a fake pair of nodes to build the exchange
                    
                    using node_t = typename required_blocks_type::node_type;
                    const node_t& src_node = src_grid.get_blocks().get_amr_node(src_lb.value);
                    const node_t& dst_node = dst_grid.get_blocks().get_amr_node(dst_lb.value);
                    
                    // Prepare for something unholy
                    node_t fake_src_node;
                    node_t fake_dst_node;
                    
                    fake_src_node.amr_position = src_node.amr_position;
                    fake_dst_node.amr_position = dst_node.amr_position;
                    
                    fake_src_node.level = src_node.level;
                    fake_dst_node.level = dst_node.level;
                    
                    fake_src_node.tag = src_node.tag;
                    fake_dst_node.tag = dst_node.tag;
                    
                    if (src_pm == 1)
                    {
                        //Add 2 to avoid the periodic neighbor conditions
                        fake_src_node.amr_position.min(match_dir).num_partitions += 2;
                        fake_src_node.amr_position.max(match_dir).num_partitions += 2;
                        
                        // "expand the grid to include the new node"
                        fake_dst_node.amr_position.min(match_dir).num_partitions = 
                            fake_src_node.amr_position.min(match_dir).num_partitions;
                    
                        fake_dst_node.amr_position.max(match_dir).num_partitions = 
                            fake_src_node.amr_position.max(match_dir).num_partitions;
                        
                        // "situate the new node at the end"
                        fake_dst_node.amr_position.min(match_dir).partition = 
                            fake_src_node.amr_position.max(match_dir).partition;
                            
                        int part_diff = dst_node.amr_position.max(match_dir).partition
                            - dst_node.amr_position.min(match_dir).partition;
                        
                        fake_dst_node.amr_position.max(match_dir).partition = 
                            fake_dst_node.amr_position.min(match_dir).partition + part_diff;
                    }
                    if (dst_pm == 1)
                    {
                        throw except::sp_exception("haven't implemented \"upwind\" interface (and this implementation is awful)");
                    }
                    
                    std::vector<node_t> proxy{fake_src_node};
                    using handle_type = typename node_t::handle_type;
                    handle_type handle(&proxy, 0);
                    const bool periodic_neighs = false;
                    fake_dst_node.create_neighbor(handle, periodic_neighs);
                    for (const auto& e:fake_dst_node.neighbors)
                    {
                        const auto transaction = get_transaction(src_arr.get_num_exchange(), src_grid, dst_grid, dst_lb.value, e);
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
        }
        
        return array_exchange_t(output0, output1);
    }
}