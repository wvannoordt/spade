#pragma once
#include <vector>

#include "print.h"
namespace cvdf::array_container
{
    template <class T> concept grid_data_container = requires(T t)
    {
        typename T::value_type;
    };
    
    template <grid_data_container container_t, typename data_t>
    void fill_container(container_t& container, const data_t& data)
    {
        //todo
    }
    
    template <typename data_t> void fill_container(std::vector<data_t>& container, const data_t& data)
    {
        for (auto& d:container) d = data;
    }
    
    template <grid_data_container container_t>
    void resize_container(container_t& container, const std::size_t& new_size)
    {
        //todo
    }
    
    template <typename data_t> void resize_container(std::vector<data_t>& container, const std::size_t& new_size)
    {
        container.resize(new_size);
    }
}