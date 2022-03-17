#pragma once

#include <type_traits>
#include <concepts>
#include "attribs.h"
#include "grid.h"

namespace cvdf::algs
{
    
    namespace detail
    {
        template <grid::multiblock_grid grid_t, grid::has_centering_type<grid::node_centered> arr_t>
        _finline_ auto get_coords(const grid_t& grid, const arr_t& arr, const int& i, const int& j, const int& k, const int& lb)
        {
            return grid.node_coords(i,j,k,lb);
        }
        
        template <grid::multiblock_grid grid_t, grid::has_centering_type<grid::cell_centered> arr_t>
        _finline_ auto get_coords(const grid_t& grid, const arr_t& arr, const int& i, const int& j, const int& k, const int& lb)
        {
            return grid.cell_coords(i,j,k,lb);
        }
        
        template <class T, typename rtype> concept xyz_ijk_callable = std::is_invocable<T, ctrs::array<rtype, 3>, int, int, int, int>::value &&
        requires(T t, const ctrs::array<rtype, 3>& x, const int& i, const int& j, const int& k, const int& lb)
        {
            {t(x,i,j,k,lb)} -> ctrs::basic_array;
        };
        
        template <class T, typename rtype> concept xyz_callable     = std::is_invocable<T, ctrs::array<rtype, 3>>::value                     &&
        requires(T t, const ctrs::array<rtype, 3>& x, const int& i, const int& j, const int& k, const int& lb)
        {
            {t(x)} -> ctrs::basic_array;
        };
        
        template <class T, typename rtype> concept ijk_callable     = std::is_invocable<T, int, int, int, int>::value                        &&
        requires(T t, const ctrs::array<rtype, 3>& x, const int& i, const int& j, const int& k, const int& lb)
        {
            {t(i,j,k,lb)} -> ctrs::basic_array;
        };
        
        template <grid::multiblock_grid grid_t, xyz_ijk_callable<typename grid_t::dtype> func_t, ctrs::vec_nd<3, double> vec_t>
        _finline_ auto forward_fillable_args(const grid_t& grid, const func_t& func, const vec_t& x, const int& i, const int& j, const int& k, const int& lb)
        {
            return func(x,i,j,k,lb);
        }
        
        template <grid::multiblock_grid grid_t, xyz_callable<typename grid_t::dtype> func_t, ctrs::vec_nd<3, double> vec_t>
        _finline_ auto forward_fillable_args(const grid_t& grid, const func_t& func, const vec_t& x, const int& i, const int& j, const int& k, const int& lb)
        {
            return func(x);
        }
        
        template <grid::multiblock_grid grid_t, ijk_callable<typename grid_t::dtype> func_t, ctrs::vec_nd<3, double> vec_t>
        _finline_ auto forward_fillable_args(const grid_t& grid, const func_t& func, const vec_t& x, const int& i, const int& j, const int& k, const int& lb)
        {
            return func(i,j,k,lb);
        }
    }
    
    template <class array_t, class callable_t>
    void fill_array(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy=grid::include_exchanges)
    {
        const auto& grid = arr.get_grid();
        auto grid_range = grid.get_range(arr.centering_type(), exchange_policy);
        for (auto maj: range(0, arr.get_major_dims().total_size()))
        {
            for (auto i: grid_range)
            {
                auto x = detail::get_coords(grid, arr, i[0], i[1], i[2], i[3]);
                auto data = detail::forward_fillable_args(grid, func, x, i[0], i[1], i[2], i[3]);
                for (int n = 0; n < data.size(); ++n)
                {
                    arr.unwrap_idx(n, i[0], i[1], i[2], i[3], maj[0]) = data[n];
                }
            }
        }
    }
}