#pragma once

#include <algorithm>

#include "grid/transactions.h"
#include "grid/get_transaction.h"
#include "grid/exchange_config.h"
#include "grid/exchange_message.h"

namespace spade::grid
{
    namespace detail
    {
        _sp_hybrid ctrs::array<int, 4> get_nijk_from_offsets(
            const std::size_t& thread_id,
            const auto& offsets,
            const auto& exchanges,
            int mine,
            int dest,
            const auto& id_func,
            const auto& send0_recv1)
        {
            constexpr int srv = send0_recv1.value;
            const auto getv = [&](const auto& v0, const auto& v1) -> const auto&
            {
                if constexpr (srv == 0) return v0;
                else return v1;
            };
            int variety = int(crnr_transaction);
            
            const std::size_t transaction_offst  = getv(offsets.send_rank_offsets, offsets.recv_rank_offsets)[dest];
            
            const auto& accums = getv(offsets.send_elem_accum, offsets.recv_elem_accum)[dest];
            const auto& offst  = getv(offsets.send_offsets,    offsets.recv_offsets)[dest];
            while(variety > int(zface_transaction))
            {
                if (accums[variety-2] > thread_id) break;
                variety--;
            }
            
            const std::size_t eff_thread_id   = thread_id - accums[variety - 1];
            const std::size_t individual_size = offsets.individual_sizes[variety - 1];
            const int         rel_exch_idx    = int(eff_thread_id/individual_size);
            const int         abs_exch_idx    = rel_exch_idx + offst[variety - 1] + transaction_offst;
            const auto&       exch_obj        = exchanges[abs_exch_idx];
            const int         local_idx       = int(eff_thread_id % individual_size);
            const auto&       bbx             = id_func(exch_obj);
            const int         nx              = bbx.size(0);
            const int         ny              = bbx.size(1);
            const int         nz              = bbx.size(2);
            
            int tmp = local_idx;
            // NOTE: local_idx = ix + iy*nx + iz*nx*ny;
            int ix = tmp % nx;
            tmp -= ix;
            tmp /= nx;
            int iy = tmp % ny;
            tmp -= iy;
            tmp /= ny;
            int iz = tmp;
            
            return ctrs::array<int, 4>(abs_exch_idx, ix, iy, iz);
        }
        
        template <typename alias_type, typename buffer_type>
        _sp_hybrid void pack_element(const alias_type& value, buffer_type& buf, const std::size_t& idx)
        {
            const std::size_t thread_pitch = alias_type::size();
            const std::size_t var_pitch    = 1;
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int v = 0; v < alias_type::size(); ++v) buf[thread_pitch*idx + var_pitch*v] = value[v];
            }
            else
            {
                buf[idx] = value;
            }
        }
        
        template <typename alias_type, typename buffer_type>
        _sp_hybrid void unpack_element(alias_type& value, const buffer_type& buf, const std::size_t& idx)
        {
            const std::size_t thread_pitch = alias_type::size();
            const std::size_t var_pitch    = 1;
            if constexpr (ctrs::basic_array<alias_type>)
            {
                for (int v = 0; v < alias_type::size(); ++v) value[v] = buf[thread_pitch*idx + var_pitch*v];
            }
            else
            {
                value = buf[idx];
            }
        }
    }
    
    template <typename array_t>
    struct arr_exchange_t
    {
        using grid_type    = typename array_t::grid_type;
        using value_type   = typename array_t::value_type;
        using device_type  = typename array_t::device_type;
        using config_type  = exchange_config_t<grid_type>;
        using message_type = exchange_message_t<value_type, device_type>;

        config_type config;
        message_type message;
        
        //TODO: return asych handle
        void exchange(array_t& array, typename grid_type::group_type& group)
        {
            using data_type  = typename array_t::value_type;
            using alias_type = typename array_t::alias_type;
            //First we need to pack all of the messages into the buffers in the message.
            auto array_img = array.image();
            const int rank = group.rank();
            
            // const auto& thing = config.send_data[0_c];
            // for (const auto& t:thing)
            // {
            //     auto lbloc  = utils::tag[partition::local](t.dest.min(3));
            //     auto lbglob = array.get_grid().get_partition().to_global(lbloc);
            //     if (lbglob.value == 6859 && t.dest.min(2) < 0)// && t.dest.min(0) == 0 && t.dest.min(1) == 0)
            //     {
            //         print(t.sender(), t.receiver());
            //         print(group.rank());
            //         print(t.dest);
            //         print(t.source);
            //     }
            // }
            // print("pause", utils::where());
            // std::cin.get();
            
            for (int dest = 0; dest < group.size(); ++dest)
            {
                bool can_short_circuit = (rank == dest);
                
                //Get the size of the injection and interpolation exchanges (in elements)
                std::size_t num_injec_cells = config.injec_offsets.send_message_size[dest];
                std::size_t num_intrp_cells = config.intrp_offsets.send_message_size[dest];
                
                // First we do all of the direct injections since they are trivial
                const auto inject_offsets = config.injec_offsets.image(array.device());
                const auto inject_exchgs  = utils::make_vec_image(config.send_data[0_c].data(array.device()));
                auto dest_buffer = utils::make_vec_image(message.send_buffers[dest].data(array.device()));
                
                
                //The first section of the buffer is only for the direct injections
                dest_buffer.csize -= alias_type::size()*num_intrp_cells;
                auto injec_range   = dispatch::ranges::make_range(0UL, num_injec_cells);
                auto injec_load    = [=] _sp_hybrid (const std::size_t& idx) mutable
                {
                    // We use this horrific function that gets the i, j, k relative to the index bounding box of the
                    // source
                    
                    const auto nijk = detail::get_nijk_from_offsets(
                        idx,
                        inject_offsets,
                        inject_exchgs,
                        rank,
                        dest,
                        [&](const auto& exch) { return exch.source; }, 0_c);
                    
                    
                    const auto& exch_obj = inject_exchgs[nijk[0]];
                    const int ix         = nijk[1];
                    const int iy         = nijk[2];
                    const int iz         = nijk[3];
                    const auto& bbx      = exch_obj.source;
                    cell_idx_t i_cell_donor
                    (
                        bbx.min(0) + ix, // i
                        bbx.min(1) + iy, // j
                        bbx.min(2) + iz, // k
                        bbx.min(3)       // lb
                    );
                    alias_type value = array_img.get_elem(i_cell_donor);
                    detail::pack_element(value, dest_buffer, idx);
                };
                
                auto injec_load_short = [=] _sp_hybrid (const std::size_t& idx) mutable
                {
                    // We use this horrific function that gets the i, j, k relative to the index bounding box of the
                    // source
                    
                    const auto nijk = detail::get_nijk_from_offsets(
                        idx,
                        inject_offsets,
                        inject_exchgs,
                        rank,
                        dest,
                        [&](const auto& exch) { return exch.source; }, 0_c);
                    
                    
                    const auto& exch_obj = inject_exchgs[nijk[0]];
                    const int ix         = nijk[1];
                    const int iy         = nijk[2];
                    const int iz         = nijk[3];
                    const auto& bbx      = exch_obj.source;
                    const auto& dstbbx   = exch_obj.dest;
                    cell_idx_t i_cell_donor
                    (
                        bbx.min(0) + ix, // i
                        bbx.min(1) + iy, // j
                        bbx.min(2) + iz, // k
                        bbx.min(3)       // lb
                    );
                    cell_idx_t i_cell_dest
                    (
                        dstbbx.min(0) + ix, // i
                        dstbbx.min(1) + iy, // j
                        dstbbx.min(2) + iz, // k
                        dstbbx.min(3)       // lb
                    );
                    array_img.set_elem(i_cell_dest, array_img.get_elem(i_cell_donor));
                };
                if (can_short_circuit) dispatch::execute(injec_range, injec_load_short, array.device());
                else                   dispatch::execute(injec_range, injec_load,       array.device());
                
                
                // Now we perform the interpolation exchanges.
                // Note that this is ruthlessly hacky.
                dest_buffer.ptr   += alias_type::size()*num_injec_cells;
                dest_buffer.csize  = alias_type::size()*num_intrp_cells;
                
                const auto interp_offsets = config.intrp_offsets.image(array.device());
                const auto interp_exchgs  = utils::make_vec_image(config.send_data[1_c].data(array.device()));
                auto intrp_range = dispatch::ranges::make_range(0UL, num_intrp_cells);
                auto intrp_load  = [=] _sp_hybrid (const std::size_t& idx) mutable
                {
                    const auto nijk = detail::get_nijk_from_offsets(
                        idx,
                        interp_offsets,
                        interp_exchgs,
                        rank,
                        dest,
                        [&](const auto& exch) { return exch.patches.dest; }, 0_c);
                    
                    const auto& exch_obj = interp_exchgs[nijk[0]];
                    const int ix         = nijk[1];
                    const int iy         = nijk[2];
                    const int iz         = nijk[3];
                    const auto& bbx      = exch_obj.patches.source;
                    constexpr int max_size = exch_obj.max_size;
                    
                    alias_type value          = data_type(0.0);
                    constexpr data_type coeff = data_type(1.0)/max_size;
                    algs::static_for<0,max_size>([&](const auto iii)
                    {
                        constexpr int iv = iii.value;
                        constexpr int d0 = (iv >> 0)&1;
                        constexpr int d1 = (iv >> 1)&1;
                        constexpr int d2 = (iv >> 2)&1;
                        int i0 = ix<<(exch_obj.i_coeff[0]+1);
                        int j0 = iy<<(exch_obj.i_coeff[1]+1);
                        int k0 = iz<<(exch_obj.i_coeff[2]+1);
                        i0 = i0 >> 1;
                        j0 = j0 >> 1;
                        k0 = k0 >> 1;
                        const int donor_di = i0 + d0*exch_obj.i_incr[0];
                        const int donor_dj = j0 + d1*exch_obj.i_incr[1];
                        const int donor_dk = k0 + d2*exch_obj.i_incr[2];
                        grid::cell_idx_t cell_idx(
                            bbx.min(0) + donor_di,
                            bbx.min(1) + donor_dj,
                            bbx.min(2) + donor_dk,
                            bbx.min(3)
                        );
                        value += array_img.get_elem(cell_idx);
                    });
                    value *= coeff;
                    detail::pack_element(value, dest_buffer, idx);
                };
                
                auto intrp_load_short  = [=] _sp_hybrid (const std::size_t& idx) mutable
                {
                    const auto nijk = detail::get_nijk_from_offsets(
                        idx,
                        interp_offsets,
                        interp_exchgs,
                        rank,
                        dest,
                        [&](const auto& exch) { return exch.patches.dest; }, 0_c);
                    
                    const auto& exch_obj = interp_exchgs[nijk[0]];
                    const int ix         = nijk[1];
                    const int iy         = nijk[2];
                    const int iz         = nijk[3];
                    const auto& bbx      = exch_obj.patches.source;
                    const auto& dst      = exch_obj.patches.dest;
                    constexpr int max_size = exch_obj.max_size;
                    
                    alias_type value          = data_type(0.0);
                    constexpr data_type coeff = data_type(1.0)/max_size;
                    algs::static_for<0,max_size>([&](const auto iii)
                    {
                        constexpr int iv = iii.value;
                        constexpr int d0 = (iv >> 0)&1;
                        constexpr int d1 = (iv >> 1)&1;
                        constexpr int d2 = (iv >> 2)&1;
                        int i0 = ix<<(exch_obj.i_coeff[0]+1);
                        int j0 = iy<<(exch_obj.i_coeff[1]+1);
                        int k0 = iz<<(exch_obj.i_coeff[2]+1);
                        i0 = i0 >> 1;
                        j0 = j0 >> 1;
                        k0 = k0 >> 1;
                        const int donor_di = i0 + d0*exch_obj.i_incr[0];
                        const int donor_dj = j0 + d1*exch_obj.i_incr[1];
                        const int donor_dk = k0 + d2*exch_obj.i_incr[2];
                        grid::cell_idx_t cell_idx(
                            bbx.min(0) + donor_di,
                            bbx.min(1) + donor_dj,
                            bbx.min(2) + donor_dk,
                            bbx.min(3)
                        );
                        value += array_img.get_elem(cell_idx);
                    });
                    value *= coeff;
                    grid::cell_idx_t idx_dest(
                            dst.min(0) + ix,
                            dst.min(1) + iy,
                            dst.min(2) + iz,
                            dst.min(3)
                        );
                    array_img.set_elem(idx_dest, value);
                };
                
                if (can_short_circuit) dispatch::execute(intrp_range, intrp_load_short, array.device());
                else                   dispatch::execute(intrp_range, intrp_load, array.device());
            }
            // At this point, the send buffers are filled.
            // Now, we handle the communication for the receive buffers.
            
            // send message
            message.send_all(group);
            
            // Now, we have to deserialize the received data. This largely looks like the serialization procedure.
            // Ideally, this can be done in its own function
            for (int dest = 0; dest < group.size(); ++dest)
            {
                bool can_short_circuit = (rank == dest);
                
                //Get the size of the injection and interpolation exchanges (in elements)
                std::size_t num_injec_cells = config.injec_offsets.recv_message_size[dest];
                std::size_t num_intrp_cells = config.intrp_offsets.recv_message_size[dest];
                
                // First we do all of the direct injections since they are trivial
                const auto inject_offsets = config.injec_offsets.image(array.device());
                const auto inject_exchgs  = utils::make_vec_image(config.recv_data[0_c].data(array.device()));
                auto src_buffer = utils::make_vec_image(message.recv_buffers[dest].data(array.device()));
                
                //The first section of the buffer is only for the direct injections
                src_buffer.csize  -= alias_type::size()*num_intrp_cells;
                auto injec_range   = dispatch::ranges::make_range(0UL, num_injec_cells);
                auto injec_load    = [=] _sp_hybrid (const std::size_t& idx) mutable
                {
                    const auto nijk = detail::get_nijk_from_offsets(
                        idx,
                        inject_offsets,
                        inject_exchgs,
                        rank,
                        dest,
                        [&](const auto& exch) { return exch.dest; }, 1_c);
                    
                    
                    const auto& exch_obj = inject_exchgs[nijk[0]];
                    const int ix         = nijk[1];
                    const int iy         = nijk[2];
                    const int iz         = nijk[3];
                    const auto& bbx      = exch_obj.dest;
                    cell_idx_t i_cell_dest
                    (
                        bbx.min(0) + ix, // i
                        bbx.min(1) + iy, // j
                        bbx.min(2) + iz, // k
                        bbx.min(3)       // lb
                    );
                    
                    alias_type value;
                    detail::unpack_element(value, src_buffer, idx);
                    array_img.set_elem(i_cell_dest, value);
                };
                
                if (!can_short_circuit) dispatch::execute(injec_range, injec_load, array.device());
                
                // Now we perform the interpolation exchanges.
                // Note that this is ruthlessly hacky.
                src_buffer.ptr   += alias_type::size()*num_injec_cells;
                src_buffer.csize  = alias_type::size()*num_intrp_cells;
                
                const auto interp_offsets = config.intrp_offsets.image(array.device());
                const auto interp_exchgs  = utils::make_vec_image(config.recv_data[1_c].data(array.device()));
                auto intrp_range = dispatch::ranges::make_range(0UL, num_intrp_cells);
                auto intrp_load  = [=] _sp_hybrid (const std::size_t& idx) mutable
                {
                        const auto nijk = detail::get_nijk_from_offsets(
                        idx,
                        interp_offsets,
                        interp_exchgs,
                        rank,
                        dest,
                        [&](const auto& exch) { return exch.patches.dest; }, 1_c);
                    
                    
                    const auto& exch_obj = interp_exchgs[nijk[0]];
                    const int ix         = nijk[1];
                    const int iy         = nijk[2];
                    const int iz         = nijk[3];
                    const auto& bbx      = exch_obj.patches.dest;
                    cell_idx_t i_cell_dest
                    (
                        bbx.min(0) + ix, // i
                        bbx.min(1) + iy, // j
                        bbx.min(2) + iz, // k
                        bbx.min(3)       // lb
                    );
                    
                    alias_type value;
                    detail::unpack_element(value, src_buffer, idx);
                    array_img.set_elem(i_cell_dest, value);
                };
                
                if (!can_short_circuit) dispatch::execute(intrp_range, intrp_load, array.device());
            }
        }
    };
    
    
    template <typename array_t>
    inline arr_exchange_t<array_t>
    make_exchange(array_t& array, ctrs::array<bool, array_t::dim()>& periodic)    
    {
        exchange_config_t<typename array_t::grid_type> config = get_exchg_config(array, periodic);
        exchange_message_t<typename array_t::value_type, typename array_t::device_type> message = get_exchg_message(array, config);
        return arr_exchange_t<array_t>{config, message};
    }
}