#pragma once

#include "sym/symbol.h"
#include "core/vec_image.h"
#include "dispatch/device_type.h"

namespace spade::ode
{
    //Represents a loose view of data stored somewhere
    template <typename data_t, sym::is_symbol mesh_sym_t, sym::is_vector variables_t>
    struct ode_buffer_t
    {
        using variable_list = variables_t;
        using mesh_symbol   = mesh_sym_t;
        
        constexpr static int num_buffers = 1 + variable_list::size();
        
        using buf_type = utils::const_vec_image_t<data_t>;
        ctrs::array<buf_type, num_buffers> buffers;
        
        template <sym::is_symbol sym_t> _sp_hybrid
        const buf_type& operator [] (const sym_t&) const
        {
            constexpr int i = 1 + variables_t::template index_of<sym_t>;
            return buffers[i];
        }
        
        _sp_hybrid
        const buf_type& operator [] (const mesh_sym_t&) const
        {
            return buffers[0];
        }
    };
    
    template <typename data_t, typename mesh_t, sym::is_vector variables_t, typename device_t>
    struct system_t
    {
        using fundamental_type = data_t;
        using value_type       = data_t;
        using mesh_type        = mesh_t;
        using variable_list    = variables_t;
        using device_type      = device_t;
        using mesh_symbol      = typename mesh_t::symbol_type;
        
        using cpu_container_type = std::vector<value_type>;
        using gpu_container_type = device::device_vector<value_type>;
        using container_type     = typename std::conditional<device::is_gpu<device_t>, gpu_container_type, cpu_container_type>::type;
        
        mesh_type mesh;
        container_type data;
        device_type device;
        
        using buffer_type = ode_buffer_t<data_t, mesh_symbol, variable_list>;
        
        system_t(const data_t&, const mesh_t& mesh_in, const variables_t&, const device_t& device_in)
        : device{device_in}, mesh{mesh_in}
        {
            const std::size_t total_elems = mesh.total_elements()*variables_t::size();
            data.resize(total_elems);
        }
    };
}