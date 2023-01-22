#pragma once

#include <type_traits>

#include "core/config.h"
#include "core/aliases.h" 
#include "core/static_for.h"

namespace spade::grid
{
    template <const int i0, const int i1> struct static_dim_t
    {
        static constexpr std::size_t size() {return i1-i0;}
        static constexpr int start()        {return i0;}
        static constexpr int end()          {return i1;}
        static constexpr int rank()         {return 1;}
        static constexpr int num_coeffs()  {return 1;}
    };

    template <typename integral_t=int> struct dynamic_dim_t
    {
        integral_t i0, i1;
        std::size_t size()  const           {return i1-i0;}
        integral_t  start() const           {return i0;}
        integral_t  end()   const           {return i1;}
        static constexpr int rank()        {return 1;}
        static constexpr int num_coeffs()  {return 1;}
        dynamic_dim_t(){}
        dynamic_dim_t(integral_t i0_in, integral_t i1_in):i0{i0_in},i1{i1_in} {}
    };

    namespace detail
    {
        template <typename dim_t, typename... dims_t> struct rank_sum_t
        {
            constexpr static int value = dim_t::rank() + rank_sum_t<dims_t...>::value;
        };
        
        template <typename dim_t> struct rank_sum_t<dim_t> {constexpr static int value = dim_t::rank();};
        
        template <typename dim_t, typename... dims_t> struct coeff_count_sum_t
        {
            constexpr static int value = dim_t::num_coeffs() + coeff_count_sum_t<dims_t...>::value;
        };
        
        template <typename dim_t> struct coeff_count_sum_t<dim_t> {constexpr static int value = dim_t::num_coeffs();};
    }

    using offset_type = long int;

    template <typename... dims_t> struct recti_view_t
    {
        constexpr static int rank()        {return detail::rank_sum_t<dims_t...>::value;}
        constexpr static int num_views()   {return sizeof...(dims_t);}
        constexpr static int num_coeffs() {return detail::coeff_count_sum_t<dims_t...>::value;}
        std::tuple<dims_t...> views;
        recti_view_t(){}
        recti_view_t(dims_t... views_in)
        {
            views = std::make_tuple(views_in...);
        }
        
        using coeff_array_t = ctrs::array<int, num_coeffs()>;
        coeff_array_t get_extent_array() const { return get_extent_array([](const auto& vw) -> auto {return vw.size();});}
        template <typename accessor_t> coeff_array_t get_extent_array(const accessor_t& accessor) const
        {
            coeff_array_t output;
            int idx = 0;
            algs::static_for<0,num_views()>([&](const auto& ii) -> void
            {
                const int i = ii.value;
                const auto& view = std::get<i>(views);
                if constexpr (requires{view.size();})
                {
                    output[idx++] = accessor(view);
                }
                else
                {
                    auto sub_coeffs = view.get_extent_array(accessor);
                    for (const auto& c: sub_coeffs) output[idx++] = c;
                }
            });
            return output;
        }
    };
    template <typename map_t> struct mem_map_t
    {
        map_t map;
        ctrs::array<std::size_t, map_t::num_coeffs()> i_coeff;
        int offset_base;
        mem_map_t(){}
        mem_map_t(const map_t& map_in):map{map_in}
        {
            compute_coeffs();
            compute_offset_base();        
        }
        
        void compute_coeffs()
        {
            auto sizes = map.get_extent_array();
            for (int i = 0; i < i_coeff.size(); ++i)
            {
                i_coeff[i] = 1;
                for (int j = 0; j < i; ++j) i_coeff[i] *= sizes[j];
            }
        }
        
        void compute_offset_base()
        {
            offset_base = 0;
            auto sizes = map.get_extent_array([](const auto& vw) -> auto {return -vw.start();});
            int idx = 0;
            for (const auto& ii: sizes) offset_base += i_coeff[idx++]*ii;
        }
        
        template <
            const int ar_idx,
            const int coeff_idx,
            typename offset_t,
            typename idx_t, typename... idxs_t>
        void rec_off_calc(offset_t& offset, const idx_t& idx, const idxs_t&... idxs) const
        {
            if constexpr (ctrs::basic_array<idx_t>)
            {
                offset += idx[ar_idx]*i_coeff[coeff_idx];
                if constexpr(ar_idx == idx_t::size()-1)
                {
                    if constexpr (sizeof...(idxs) > 0)
                    {
                        rec_off_calc<0,coeff_idx+1>(offset, idxs...);
                    }
                }
                else
                {
                    rec_off_calc<ar_idx+1,coeff_idx+1>(offset, idx, idxs...);
                }
            }
            else
            {
                offset += idx*i_coeff[coeff_idx];
                if constexpr (sizeof...(idxs) > 0)
                {
                    rec_off_calc<0,coeff_idx+1>(offset, idxs...);
                }
            }
        }
        
        template <typename... idxs_t>
        requires (true)
        offset_type compute_offset(const idxs_t&... idxs) const
        {
            offset_type output = 0;
            rec_off_calc<0,0>(output, idxs...);
            return output + offset_base;
        }
    };
}