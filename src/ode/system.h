#pragma once

#include "sym/symbol.h"
#include "core/vec_image.h"
#include "dispatch/device_type.h"

namespace spade::ode
{
    //Note: this implementation is pretty much garbage
    
    //Represents a loose view of data stored somewhere
    template <typename container_img_t, typename data_t, sym::is_symbol mesh_sym_t, sym::is_vector variables_t>
    struct ode_buffer_t
    {
        using variable_list = variables_t;
        using mesh_symbol   = mesh_sym_t;
        
        constexpr static int num_buffers = 1 + variable_list::size();
        
        using buf_type = container_img_t;
        ctrs::array<buf_type, num_buffers> buffers;
        
        template <sym::is_symbol sym_t> _sp_hybrid
        const buf_type& operator [] (const sym_t&) const
        {
            constexpr int i = 1 + variables_t::template index_of<sym_t>;
            return buffers[i];
        }
        
        template <sym::is_symbol sym_t> _sp_hybrid
        buf_type& operator [] (const sym_t&)
        {
            constexpr int i = 1 + variables_t::template index_of<sym_t>;
            return buffers[i];
        }
        
        _sp_hybrid const buf_type& operator [] (const mesh_sym_t&) const
        {
            return buffers[0];
        }
        
        _sp_hybrid buf_type& operator [] (const mesh_sym_t&)
        {
            return buffers[0];
        }
    };
    
    
    // Todo: improve this, this is a little bit janky
    template <typename container_image_t, typename data_t, typename mesh_t, sym::is_vector variables_t>
    struct system_image_t
    {
        container_image_t ys_raw;
        container_image_t data_raw;
        using buffer_type       = ode_buffer_t<utils::vec_image_t<data_t>,       data_t, typename mesh_t::symbol_type, variables_t>;
        using const_buffer_type = ode_buffer_t<utils::const_vec_image_t<data_t>, data_t, typename mesh_t::symbol_type, variables_t>;
        
        _sp_hybrid buffer_type get_buffer(const std::size_t& instance)
        {
            buffer_type output;
            const std::size_t ngrd = ys_raw.size();
            
            output.buffers[0] = container_image_t{ys_raw.ptr + ngrd*instance, ngrd};
            for (int j = 0; j < output.buffers.size() - 1; ++j)
            {
                output.buffers[1+j] = container_image_t{data_raw.ptr + ngrd*instance*j, ngrd};
            }
            
            return output;
        }
        
        _sp_hybrid const_buffer_type get_buffer(const std::size_t& instance) const
        {
            const_buffer_type output;
            const std::size_t ngrd = ys_raw.size();
            
            output.buffers[0] = container_image_t{ys_raw.ptr + ngrd*instance, ngrd};
            for (int j = 0; j < output.buffers.size() - 1; ++j)
            {
                output.buffers[1+j] = container_image_t{data_raw.ptr + ngrd*instance*j, ngrd};
            }
            
            return output;
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
        using container_type     = device::auto_vector<value_type, device_t>;
        
        mesh_type mesh;
        container_type data;
        device_type devc;
        
        using buffer_type       = ode_buffer_t<utils::vec_image_t<data_t>,       data_t, typename mesh_t::symbol_type, variables_t>;
        using const_buffer_type = ode_buffer_t<utils::const_vec_image_t<data_t>, data_t, typename mesh_t::symbol_type, variables_t>;
        
        system_t(const data_t&, const mesh_t& mesh_in, const variables_t&, const device_t& device_in)
        : devc{device_in}, mesh{mesh_in}
        {
            const std::size_t total_elems = mesh.total_elements()*variables_t::size();
            data.resize(total_elems);
        }
        
        device_type device() const { return devc; }
        
        using image_type       = system_image_t<utils::vec_image_t<data_t>, data_t, mesh_t, variables_t>;
        using const_image_type = system_image_t<utils::const_vec_image_t<data_t>, data_t, mesh_t, variables_t>;
        image_type image()
        {
            return image_type{utils::make_vec_image(mesh.xs), utils::make_vec_image(data)};
        }
        const const_image_type image() const
        {
            return const_image_type{utils::make_vec_image(mesh.xs), utils::make_vec_image(data)};
        }
    };
}