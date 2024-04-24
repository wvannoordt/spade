#pragma once

#include "grid/exchange_config.h"

namespace spade::grid
{
    template <typename data_t, typename device_t>
    struct exchange_message_t
    {
        std::vector<device::shared_vector<data_t, device::pinned_allocator_t<data_t>, device::device_allocator_t<data_t>>> send_buffers;
        std::vector<device::shared_vector<data_t, device::pinned_allocator_t<data_t>, device::device_allocator_t<data_t>>> recv_buffers;
        
        template <typename group_t>
        void send_all(const group_t& group)
        {
            //Do some logic here!
            const auto dev = device_t();
            // const auto dev = device::cpu;
            
            if constexpr (device::is_gpu<device_t>)
            {
                int rank = 0;
                for (auto& s:send_buffers)
                {
                    const auto channel = parallel::channel_from_device(dev);
                    if (channel == parallel::cpu_messg) s.itransfer();
                    ++rank;
                }
            }
            
            using message_type = utils::vec_image_t<data_t>;
            for (int p = 0; p < group.size(); ++p)
            {
                const auto channel = parallel::channel_from_device(dev);
                // const auto dev = device::cpu;
                auto send_buf = utils::make_vec_image(send_buffers[p].data(dev));
                auto recv_buf = utils::make_vec_image(recv_buffers[p].data(dev));
                group.post_send(send_buf, p, channel);
                group.post_recv(recv_buf, p);
            }
            
            group.template send_all<message_type>();
            
            if constexpr (device::is_gpu<device_t>)
            {
                int rank = 0;
                for (auto& r:recv_buffers)
                {
                    const auto channel = parallel::channel_from_device(dev);
                    if (channel == parallel::cpu_messg) r.transfer();
                    ++rank;
                }
            }
        }
        
        // Used to resize receive buffers to the correct size 
        template <typename group_t>
        void assure_recv_buf_size(const group_t& group)
        {
            group.sync();
            using tmp_dev_t = device::cpu_t;
            exchange_message_t<std::size_t, tmp_dev_t> size_msg;
            size_msg.set_size(group.size());
            for (int i = 0; i < group.size(); ++i)
            {
                size_msg.send_buffers[i].resize(1);
                size_msg.send_buffers[i][0] = send_buffers[i].size();
                size_msg.recv_buffers[i].resize(1);
            }
            size_msg.send_all(group);
            group.sync();
            for (int i = 0; i < group.size(); ++i)
            {
                std::size_t req_size = size_msg.recv_buffers[i][0];
                recv_buffers[i].resize(req_size);
            }
            group.sync();
        }
        
        void set_size(const std::size_t& pool_size)
        {
            send_buffers.resize(pool_size);
            recv_buffers.resize(pool_size);
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
        output.set_size(group.size());
        
        for (int prc = 0; prc < group.size(); ++prc)
        {
            if (prc != group.rank())
            {
                output.send_buffers[prc].resize(num_vars*(g_config.injec_offsets.send_message_size[prc] + g_config.intrp_offsets.send_message_size[prc]));
                output.recv_buffers[prc].resize(num_vars*(g_config.injec_offsets.recv_message_size[prc] + g_config.intrp_offsets.recv_message_size[prc]));
                output.send_buffers[prc].transfer();
                output.recv_buffers[prc].transfer();
            }
        }
        
        return output;
    }
}