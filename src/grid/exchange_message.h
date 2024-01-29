#pragma once

#include "grid/exchange_config.h"

namespace spade::grid
{
    template <typename data_t, typename device_t>
    struct exchange_message_t
    {
        std::vector<device::shared_vector<data_t>> send_buffers;
        std::vector<device::shared_vector<data_t>> recv_buffers;
        
        template <typename group_t>
        void send_all(group_t& group)
        {
            //Horribly inefficient implementation!
            for (auto& s:send_buffers) s.itransfer();
            
            //Let's forget about optimizations for now!
            using message_type = utils::vec_image_t<data_t>;
            for (int p = 0; p < group.size(); ++p)
            {
                auto send_buf = utils::make_vec_image(send_buffers[p].data(device::cpu));
                auto recv_buf = utils::make_vec_image(recv_buffers[p].data(device::cpu));
                group.post_send(send_buf, p);
                group.post_recv(recv_buf, p);
            }
            
            group.template send_all<message_type>();
            
            for (auto& r:recv_buffers) r.transfer();
        }
    };
    
    template <typename array_t, typename config_t>
    inline exchange_message_t<typename array_t::value_type, typename array_t::device_type>
    get_exchg_message(array_t& array, const config_t& g_config)
    {
        auto& group = array.get_grid().group();
        constexpr int num_vars = array_t::alias_type::size();
        using output_t = exchange_message_t<typename array_t::value_type, typename array_t::device_type>;
        output_t output;
        
        output.send_buffers.resize(group.size());
        output.recv_buffers.resize(group.size());
        
        for (int prc = 0; prc < group.size(); ++prc)
        {
            output.send_buffers[prc].resize(num_vars*(g_config.injec_offsets.send_message_size[prc] + g_config.intrp_offsets.send_message_size[prc]));
            output.recv_buffers[prc].resize(num_vars*(g_config.injec_offsets.recv_message_size[prc] + g_config.intrp_offsets.recv_message_size[prc]));
            output.send_buffers[prc].transfer();
            output.recv_buffers[prc].transfer();
        }
        
        return output;
    }
}