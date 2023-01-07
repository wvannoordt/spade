#pragma once

#include <iostream>

namespace spade::array_containers
{
    template <typename T> concept is_array_container = requires(T t, const typename T::value_type& v, const typename T::size_type& n)
    {
        t.resize(n);
        t.resize(n, v);
        t.size();
        t.fill(v);
        t[n];
    };
    
    template <typename data_t, typename derived_t> struct array_container
    {
        using size_type = std::size_t;
        using value_type = data_t;
        using self_type = derived_t;
        self_type& self() { return *(static_cast<self_type*>(this)); }
        
        //required interface implementations
        virtual void resize(const size_type& new_size, const value_type val = data_t()) = 0;
        virtual size_type size() const = 0;
        virtual void fill(const data_t& val) = 0;
        virtual value_type& operator [] (const size_type i) = 0;
        virtual const value_type& operator [] (const size_type i) const = 0;
        
        //need something about capture types here (TBD)
    };
}