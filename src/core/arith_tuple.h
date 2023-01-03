#pragma once

#include <utility>

#include "core/aliases.h"
#include "core/static_for.h"
#include "core/get.h"

namespace spade::ctrs
{
    namespace detail
    {
        //This is a mess...
        template <typename T> concept has_data = requires(const T& t) {t.data;};
        
        template <const int i, typename rhs_t> const auto& access_const_elem(const rhs_t& rhs) {return rhs;}
        
        template <const int i, typename... elems_t>
        const auto& access_const_elem(const aliases::tuple<elems_t...>& rhs) {return get<i>(rhs);}
        
        template <const int i, has_data rhs_t>
        const auto& access_const_elem(const rhs_t& rhs) {return get<i>(rhs.data);}
    }
    
    template <typename... elems_t> struct arith_tuple
    {
        aliases::tuple<elems_t...> data;
        arith_tuple(const elems_t&... elems)
        {
            data = std::make_tuple(elems...);
        }
        
        arith_tuple(elems_t&&... elems)
        {
            data = std::make_tuple(elems...);
        }
        
        auto& get_tuple() {return data;}
        const auto& get_tuple() const {return data;}
        
        template <typename rhs_t>
        arith_tuple& operator += (const rhs_t& rhs)
        {
            algs::static_for<0,sizeof...(elems_t)>([&](const auto& i) -> void
            {
                const int ii = i.value;
                get<ii>(data) += detail::access_const_elem<ii>(rhs);
            });
            return *this;
        }
        
        template <typename rhs_t>
        arith_tuple& operator -= (const rhs_t& rhs)
        {
            algs::static_for<0,sizeof...(elems_t)>([&](const auto& i) -> void
            {
                const int ii = i.value;
                get<ii>(data) -= detail::access_const_elem<ii>(rhs);
            });
            return *this;
        }
        
        template <typename rhs_t>
        arith_tuple& operator *= (const rhs_t& rhs)
        {
            algs::static_for<0,sizeof...(elems_t)>([&](const auto& i) -> void
            {
                const int ii = i.value;
                get<ii>(data) *= detail::access_const_elem<ii>(rhs);
            });
            return *this;
        }
        
        template <typename rhs_t>
        arith_tuple& operator /= (const rhs_t& rhs)
        {
            algs::static_for<0,sizeof...(elems_t)>([&](const auto& i) -> void
            {
                const int ii = i.value;
                get<ii>(data) *= detail::access_const_elem<ii>(rhs);
            });
            return *this;
        }
    };
}