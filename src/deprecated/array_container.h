#pragma once
#include <vector>

#include "core/print.h"
#include "core/ctrs.h"
namespace spade::array_container
{
    template <class T> concept grid_data_container = requires(T t, std::size_t i)
    {
        t[i];
        t.size();
        typename T::value_type;
    };
    
    //Container filling
    template <grid_data_container container_t, typename data_t>
    static void fill_container(container_t& container, const data_t& data)
    {
        //todo
    }
    
    template <typename data_t> static void fill_container(std::vector<data_t>& container, const data_t& data)
    {
        for (auto& d:container) d = data;
    }
    
    template <typename data_t> static void fill_container(ctrs::unsafe_vector_alias_t<data_t>& container, const data_t& data)
    {
        for (auto i: range(0, container.size())) container[i] = data;
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
    
    template <typename data_t> static void resize_container(ctrs::unsafe_vector_alias_t<data_t>& container, const std::size_t& new_size)
    {
        container.safe_size = new_size;
    }
}