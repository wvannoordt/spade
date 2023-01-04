#pragma once

#include <type_traits>
#include <ratio>

#include "core/utils.h"

namespace spade::time_integration
{
    namespace detail
    {
        template <typename... values_t> struct sr_vec_t
        {
            template <const int i> using elem_t = utils::get_pack_type<i, values_t...>::type;
            static constexpr std::size_t length() {return sizeof...(values_t);}
        };
        template <typename... vecs_t> struct sr_mat_t
        {
            template <const int i> using elem_t = utils::get_pack_type<i, vecs_t...>::type;
            static constexpr std::size_t rows() {return sizeof...(vecs_t);}
            static constexpr std::size_t cols() {return elem_t<0>::length();}
        };
    }
    
    template <typename table_t, typename accum_t, typename dt_t>
    requires(table_t::rows() == dt_t::length()) // need other constraints here
    struct rk_t
    {
        using table_type = table_t; //how do we build the solution from residuals at each substep?
        using accum_type = accum_t; //How do we combine the substep residuals?
        using dt_type    = dt_t;    //At what points in time do we evaluate t?
        
        static constexpr std::size_t var_size() {return 1;}
        static constexpr std::size_t rhs_size() {return table_t::cols();}
    };
    
    using rk2_t = rk_t<
        detail::sr_mat_t<
            detail::sr_vec_t<std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1,2>, std::ratio<0>>
        >,
        detail::sr_vec_t<std::ratio<0>, std::ratio<1>>,
        detail::sr_vec_t<std::ratio<0>, std::ratio<1,2>>
    >;
    
    using rk4_t = rk_t<
        detail::sr_mat_t<
            detail::sr_vec_t<std::ratio<0>,   std::ratio<0>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1,2>, std::ratio<0>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<0>,   std::ratio<1,2>, std::ratio<1,2>, std::ratio<0>>,
            detail::sr_vec_t<std::ratio<0>,   std::ratio<0>,   std::ratio<1>,   std::ratio<0>>
        >,
        detail::sr_vec_t<std::ratio<1,6>, std::ratio<1,3>, std::ratio<1,3>, std::ratio<1,6>>,
        detail::sr_vec_t<std::ratio<0>,   std::ratio<1,2>, std::ratio<1,2>, std::ratio<1>>
    >;
    
    using ssprk34_t = rk_t<
        detail::sr_mat_t<
            detail::sr_vec_t<std::ratio<0>,   std::ratio<0>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1,2>, std::ratio<0>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1,2>, std::ratio<1,2>, std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1,6>, std::ratio<1,6>, std::ratio<1,6>, std::ratio<0>>
        >,
        detail::sr_vec_t<std::ratio<1,6>, std::ratio<1,6>, std::ratio<1,6>, std::ratio<1,2>>,
        detail::sr_vec_t<std::ratio<0>,   std::ratio<1,2>, std::ratio<1>,   std::ratio<1,2>>
    >;
    
    using rk_38_t = rk_t<
        detail::sr_mat_t<
            detail::sr_vec_t<std::ratio<0>,    std::ratio<0>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1,3>,  std::ratio<0>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<-1,3>, std::ratio<1>,   std::ratio<0>,   std::ratio<0>>,
            detail::sr_vec_t<std::ratio<1>,    std::ratio<-1>,  std::ratio<1>,   std::ratio<0>>
        >,
        detail::sr_vec_t<std::ratio<1,8>, std::ratio<3,8>, std::ratio<3,8>, std::ratio<1,8>>,
        detail::sr_vec_t<std::ratio<0>,   std::ratio<1,3>, std::ratio<2,3>, std::ratio<1>>
    >;
}