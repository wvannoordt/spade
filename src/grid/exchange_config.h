#pragma once

#include "core/parallel.h"
#include "dispatch/shared_vector.h"

namespace spade::grid
{
    template <const int variety_size, template <typename> typename container_t = device::basic_shared_vector>
    struct message_offset_t
    {
        using image_type = message_offset_t<variety_size, utils::const_vec_image_t>;
        
        container_t<std::size_t> send_rank_offsets;
        container_t<std::size_t> send_rank_sizes;
        
        container_t<std::size_t> recv_rank_offsets;
        container_t<std::size_t> recv_rank_sizes;
        
        container_t<ctrs::array<std::size_t, variety_size>> send_sizes;
        container_t<ctrs::array<std::size_t, variety_size>> send_offsets;
        
        container_t<ctrs::array<std::size_t, variety_size>> send_elem_size;
        container_t<ctrs::array<std::size_t, variety_size>> send_elem_accum;
        
        container_t<ctrs::array<std::size_t, variety_size>> recv_sizes;
        container_t<ctrs::array<std::size_t, variety_size>> recv_offsets;
        
        container_t<ctrs::array<std::size_t, variety_size>> recv_elem_size;
        container_t<ctrs::array<std::size_t, variety_size>> recv_elem_accum;
        
        ctrs::array<std::size_t, variety_size> individual_sizes;
        
        container_t<std::size_t> send_message_size;
        container_t<std::size_t> recv_message_size;
        
        template <typename device_t>
        image_type image(const device_t& device) const
        {
            // Where is static reflection??????
            image_type output;
            output.send_rank_offsets = utils::make_vec_image(send_rank_offsets.data(device));
            output.send_rank_sizes   = utils::make_vec_image(send_rank_sizes.data(device));
            
            output.recv_rank_offsets = utils::make_vec_image(recv_rank_offsets.data(device));
            output.recv_rank_sizes   = utils::make_vec_image(recv_rank_sizes.data(device));
            
            output.send_sizes     = utils::make_vec_image(send_sizes.data(device));
            output.send_offsets   = utils::make_vec_image(send_offsets.data(device));
            
            output.recv_sizes     = utils::make_vec_image(recv_sizes.data(device));
            output.recv_offsets   = utils::make_vec_image(recv_offsets.data(device));
            
            output.send_message_size = utils::make_vec_image(send_message_size.data(device));
            output.recv_message_size = utils::make_vec_image(recv_message_size.data(device));
            
            output.recv_elem_size  = utils::make_vec_image(recv_elem_size.data(device));
            output.recv_elem_accum = utils::make_vec_image(recv_elem_accum.data(device));
            
            output.send_elem_size  = utils::make_vec_image(send_elem_size.data(device));
            output.send_elem_accum = utils::make_vec_image(send_elem_accum.data(device));
            
            output.individual_sizes = individual_sizes;
            
            return output;
        }
        
        void count(const auto& group, const auto& sends, const auto& recvs)
        {
            // Note that all sends have rank as the sender and all recvs have rank as the 
            // receiver
            
            int my_rank = group.rank();
            int gp_size = group.size();
            
            const auto calc_data = [&](auto& rank_size_vec, auto& rank_ofst_vec, const auto& id_func, const auto& idata)
            {
                int last_recv = 0;
                rank_size_vec.resize(gp_size, 0);
                for (const auto& s: idata)
                {
                    int rank = id_func(s);
                    rank_size_vec[rank]++;
                }
                
                std::size_t offset_loc = 0;
                for (const auto& sz: rank_size_vec)
                {
                    rank_ofst_vec.push_back(offset_loc);
                    offset_loc += sz;
                }
            };

            calc_data(send_rank_sizes, send_rank_offsets, [](const auto& tr) { return tr.receiver(); }, sends);
            calc_data(recv_rank_sizes, recv_rank_offsets, [](const auto& tr) { return tr.sender();   }, recvs);
            
            send_sizes.resize(gp_size, 0UL);
            send_offsets.resize(gp_size, 0UL);
            send_elem_accum.resize(gp_size, 0UL);
            send_elem_size.resize(gp_size, 0UL);
            recv_sizes.resize(gp_size, 0UL);
            recv_offsets.resize(gp_size, 0UL);
            recv_elem_accum.resize(gp_size, 0UL);
            recv_elem_size.resize(gp_size, 0UL);
            
            const auto calc_variety_data = [&](
                auto& v_sizes,
                auto& v_offst,
                const auto& irank_size,
                const auto& irank_ofst,
                const auto& list)
            {
                for (int prc = 0; prc < gp_size; ++prc)
                {
                    std::size_t exch_start = irank_ofst[prc];
                    std::size_t exch_end   = irank_ofst[prc] + irank_size[prc];
                    auto& sizes_here = v_sizes[prc];
                    sizes_here = 0UL;
                    int cur_priority_level = int(crnr_transaction);
                    std::size_t cur_size = 0;
                    for (std::size_t i_exchg = exch_start; i_exchg < exch_end; ++i_exchg)
                    {
                        const auto& exch = list[i_exchg];
                        int tagval = int(exch.get_tag());
                        if (cur_priority_level != exch.get_tag())
                        {
                            sizes_here[cur_priority_level - 1] = cur_size;
                            cur_priority_level = int(exch.get_tag());
                            cur_size = 0;
                        }
                        cur_size++;
                    }
                    sizes_here[int(zface_transaction)-1] = cur_size;
                    std::size_t cur_offst = 0;
                    auto& offst_here = v_offst[prc];
                    offst_here = 0UL;
                    for (int i = 0; i < variety_size; ++i)
                    {
                        int idx = int(crnr_transaction) - i - 1;
                        offst_here[idx] = cur_offst;
                        cur_offst += sizes_here[idx];
                    }
                    
                }
            };
            calc_variety_data(send_sizes, send_offsets, send_rank_sizes, send_rank_offsets, sends);
            calc_variety_data(recv_sizes, recv_offsets, recv_rank_sizes, recv_rank_offsets, recvs);
            
            const auto calc_accums = [&](const auto& exchs, auto& elems, auto& accums, const auto& id_func)
            {
                for (const auto& exch: exchs)
                {
                    int id = id_func(exch);
                    elems[id][int(exch.get_tag())-1] += exch.recv_volume();
                    individual_sizes[int(exch.get_tag())-1] = exch.recv_volume();
                }
                for (int prc = 0; prc < elems.size(); ++prc)
                {
                    const auto& elems_loc  = elems [prc];
                    auto& accums_loc = accums[prc];
                    int vari = int(crnr_transaction);
                    std::size_t size_loc = 0;
                    for (int i = 0; i < elems_loc.size(); ++i)
                    {
                        accums_loc[vari-1] = size_loc;
                        size_loc += elems_loc[vari-1];
                        vari--;
                    }
                }
            };
            individual_sizes = 1;
            calc_accums(sends, send_elem_size, send_elem_accum, [&](const auto& thing) { return thing.receiver(); });
            calc_accums(recvs, recv_elem_size, recv_elem_accum, [&](const auto& thing) { return thing.sender(); });
            
            send_message_size.resize(gp_size, 0UL);
            recv_message_size.resize(gp_size, 0UL);
            for (const auto& recv_exch: recvs)
            {
                recv_message_size[recv_exch.sender()]   += recv_exch.recv_volume();
            }
            for (const auto& send_exch: sends)
            {
                send_message_size[send_exch.receiver()] += send_exch.recv_volume();
            }
            send_rank_offsets.transfer();
            send_rank_sizes.transfer();
            recv_rank_offsets.transfer();
            recv_rank_sizes.transfer();
            send_sizes.transfer();
            send_offsets.transfer();
            recv_sizes.transfer();
            recv_offsets.transfer();
            send_message_size.transfer();
            recv_message_size.transfer();
            send_elem_accum.transfer();
            send_elem_size.transfer();
            recv_elem_accum.transfer();
            recv_elem_size.transfer();
        }
        
        //For debugging purposes:
        std::ostream& report(std::ostream& os, const auto& sends, const auto& recvs) const
        {
            const auto dbg = [&](const auto& list, const auto& idx_i)
            {
                constexpr int sendrecv = idx_i.value;
                
                const auto elm = [&](const auto& v0, const auto& v1) -> const auto&
                {
                    if constexpr (sendrecv == 0) { return v0; }
                    else return v1;
                };
                
                const auto st  = [&](const auto& e) { return utils::pad_str(e, 6); };
                const auto st2 = [&](const auto& e) { return utils::pad_str(e, 11); };
                int idx = 0;
                for (const auto& exch: list)
                {
                    os << "ID: "        << st(idx);
                    os << "  To: "      << st(exch.receiver());
                    os << "  From: "    << st(exch.sender());
                    os << "  Size: "    << st(exch.recv_volume());
                    os << "  Variety: " << trs2str(exch.get_tag());
                    os << "\n";
                    ++idx;
                }
                int wrld_size = send_rank_offsets.size();
                os << "=================\n";
                os << "STATS:\n";
                os << "----------------\n";
                os << "Recv size by variety:\n";
                for (int i = int(crnr_transaction); i > int(null_transaction); --i)
                {
                    os << " -- " << st2(individual_sizes[i-1]) << " elems per ";
                    os << trs2str(transaction_tag_t(i)) << "\n";
                }
                os << "Offsets:\n";
                
                const auto& sizes = elm(send_rank_sizes,   recv_rank_sizes);
                const auto& offst = elm(send_rank_offsets, recv_rank_offsets);
                
                for (int prc = 0; prc < wrld_size; ++prc)
                {
                    os << "Rank: "    << st(prc);
                    os << " Size: "   << st(sizes[prc]);
                    os << " Offset: " << st(offst[prc]);
                    os << " Start: "  << st(offst[prc]);
                    os << " End: "    << st(offst[prc] + sizes[prc]) << "\n";
                    
                    const auto& vari_sizes = elm(send_sizes,   recv_sizes)[prc];
                    const auto& vari_offts = elm(send_offsets, recv_offsets)[prc];
                    
                    int vari = int(crnr_transaction);
                    os << " -- Partial: \n";
                    for (int i = 0; i < vari_sizes.size(); ++i)
                    {
                        os << " -- -- Size: " << st(vari_sizes[vari-1]);
                        os << " Offset: " << st(vari_offts[vari-1]);
                        os << " Variety: " << trs2str(transaction_tag_t(vari));
                        vari--;
                        os << "\n";
                    }
                    
                    const auto& elem_szs = elm(send_elem_size, recv_elem_size)[prc];
                    const auto& elem_acm = elm(send_elem_accum, recv_elem_accum)[prc];
                    vari = int(crnr_transaction);
                    os << " -- Elements: partial / cumulative / variety:\n";
                    for (int i = 0; i < elem_acm.size(); ++i)
                    {
                        os << " --  -- " << st(elem_szs[vari-1]);
                        os << " / " << st(elem_acm[vari-1]);
                        os << " / " << trs2str(transaction_tag_t(vari)) << "\n";
                        vari--;
                    }
                }
            };
            
            os << "Raw data: sends\n";
            dbg(sends, 0_c);
            os << "Raw data: recvs\n";
            dbg(recvs, 1_c);
            return os;
        }
    };
    
    template <typename grid_t>
    struct exchange_config_t
    {
        utils::poly_vec_t<grid_rect_copy_t, patch_fill_t<grid_t::dim()>> send_data;
        utils::poly_vec_t<grid_rect_copy_t, patch_fill_t<grid_t::dim()>> recv_data;
        int my_rank;
        
        
        constexpr static int injec_table_size = int(crnr_transaction); //See transaction_tag_t, transactions.h
        constexpr static int intrp_table_size = int(crnr_transaction); //See transaction_tag_t, transactions.h
        
        message_offset_t<injec_table_size> injec_offsets;
        message_offset_t<intrp_table_size> intrp_offsets;
        
        template <typename transaction_t>
        void add_send(const transaction_t& elem)
        {
            if (elem.sender() == my_rank)
            {
                send_data.add(elem);
            }
        }
        
        template <typename transaction_t>
        void add_recv(const transaction_t& elem)
        {
            if (elem.receiver() == my_rank)
            {
                recv_data.add(elem);
            }
        }
    };
    
    namespace detail
    {
        template <typename array_t>
        inline exchange_config_t<typename array_t::grid_type> get_exchange_config(const array_t& array, ctrs::array<bool, array_t::dim()>& periodic)
        {
            const auto& grid  = array.get_grid();
            const auto& group = grid.group();
            const auto  num_exchanges = array.get_num_exchange();
            using grid_type = typename array_t::grid_type;
            exchange_config_t<grid_type> output;
            output.my_rank = group.rank();
            
            std::size_t ct = 0;
            for (auto lb = utils::tag[partition::global](0); lb.value < grid.get_num_global_blocks(); ++lb.value)
            {
                const auto& neighs = grid.get_blocks().get_neighs(lb.value);
                // if (global::debug) mf << ct << ": " << lb.value << " --> " << neighs.size() << "\n";
                for (const auto& e:neighs)
                {                    
                    bool ignore_from_periodic = false;
                    const auto& is_dom_bdy = grid.is_domain_boundary(lb);
                    for (int d = 0; d < array_t::dim(); ++d)
                    {
                        bool loc = !periodic[d] && ((is_dom_bdy.min(d) && e.edge[d] == -1) || (is_dom_bdy.max(d) && e.edge[d] == 1));
                        ignore_from_periodic = ignore_from_periodic || loc;
                    }
                    
                    if (!ignore_from_periodic)
                    {
                        const auto transaction = get_transaction(num_exchanges, grid, grid, lb.value, e);
                        
                        // Note that we add as both and the internal logic inside
                        // these calls will handle the rank checking
                        
                        if (transaction.reducible())
                        {
                            output.add_send(transaction.reduce());
                            output.add_recv(transaction.reduce());
                        }
                        else
                        {
                            output.add_send(transaction);
                            output.add_recv(transaction);
                        }
                    }
                }
            }
            
            output.send_data[0_c].transfer();
            output.send_data[1_c].transfer();
            output.recv_data[0_c].transfer();
            output.recv_data[1_c].transfer();
            return output;
        }
    }
    
    template <typename array_t>
    inline auto get_exchg_config(const array_t& array, ctrs::array<bool, array_t::dim()>& periodic)    
    {
        const auto& group = array.get_grid().group();
        auto config = detail::get_exchange_config(array, periodic);
        
        const auto send_direct_inj_sort = [&](const auto& ex1, const auto& ex2)
        {
            if (ex1.receiver() != ex2.receiver())  return ex1.receiver() < ex2.receiver();
            if (ex1.get_tag()  != ex2.get_tag())   return ex1.get_tag()  > ex2.get_tag();
            return ex1.recv_volume()  < ex2.recv_volume(); //We should never actually get here!
        };
        
        const auto recv_direct_inj_sort = [&](const auto& ex1, const auto& ex2)
        {
            if (ex1.sender()   != ex2.sender())    return ex1.sender()   < ex2.sender();
            if (ex1.get_tag()  != ex2.get_tag())   return ex1.get_tag()  > ex2.get_tag();
            return ex1.recv_volume()  < ex2.recv_volume(); //We should never actually get here!
        };
        
        //Note that use of stable_sort is critical!
        std::stable_sort(config.send_data[0_c].begin(), config.send_data[0_c].end(), send_direct_inj_sort);
        std::stable_sort(config.recv_data[0_c].begin(), config.recv_data[0_c].end(), recv_direct_inj_sort);
        std::stable_sort(config.send_data[1_c].begin(), config.send_data[1_c].end(), send_direct_inj_sort);
        std::stable_sort(config.recv_data[1_c].begin(), config.recv_data[1_c].end(), recv_direct_inj_sort);
        
        config.send_data[0_c].transfer();
        config.recv_data[0_c].transfer();
        config.send_data[1_c].transfer();
        config.recv_data[1_c].transfer();
               
        // This appears to work now!
        config.injec_offsets.count(group, config.send_data[0_c], config.recv_data[0_c]);
        config.intrp_offsets.count(group, config.send_data[1_c], config.recv_data[1_c]);
        
        
        // std::string fname0 = utils::zfill(group.rank(), 2) + ".inject.log";
        // std::ofstream mf0(fname0);
        // config.injec_offsets.report(mf0, config.send_data[0_c], config.recv_data[0_c]);
        
        // std::string fname1 = utils::zfill(group.rank(), 2) + ".interp.log";
        // std::ofstream mf1(fname1);
        // config.intrp_offsets.report(mf1, config.send_data[1_c], config.recv_data[1_c]);
        
        group.sync();
        return config;
    }
}