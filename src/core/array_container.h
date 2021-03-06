#pragma once
#include <vector>

#include "core/print.h"
namespace spade::array_container
{
    template <class T> concept grid_data_container = requires(T t, std::size_t i)
    {
        t[i];
        t.size();
        typename T::value_type;
    };
    
    template <grid_data_container container_t, typename data_t>
    static void fill_container(container_t& container, const data_t& data)
    {
        //todo
    }
    
    template <typename data_t> static void fill_container(std::vector<data_t>& container, const data_t& data)
    {
        for (auto& d:container) d = data;
    }
    
    template <grid_data_container container_t>
    void resize_container(container_t& container, const std::size_t& new_size)
    {
        //todo
    }
    
    template <typename data_t> static void resize_container(std::vector<data_t>& container, const std::size_t& new_size)
    {
        container.resize(new_size);
    }
}