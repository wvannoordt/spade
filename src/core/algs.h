#pragma once

#include <type_traits>
#include <concepts>

#include "core/attribs.h"
#include "core/grid.h"
#include "core/reduce_ops.h"
#include "core/md_loop.h"

namespace spade::algs
{
    
    namespace detail
    {
        template <grid::multiblock_grid grid_t, grid::has_centering_type<grid::node_centered> arr_t>
        _finline_ auto get_coords(const grid_t& grid, const arr_t& arr, const int& i, const int& j, const int& k, const int& lb)
        {
            grid::node_idx_t i_n(i, j, k, lb);
            return grid.get_coords(i_n);
        }
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <grid::multiblock_grid grid_t, grid::has_centering_type<grid::cell_centered> arr_t>
        _finline_ auto get_coords(const grid_t& grid, const arr_t& arr, const int& i, const int& j, const int& k, const int& lb)
        {
            grid::cell_idx_t i_c(i, j, k, lb);
            return grid.get_coords(i_c);
        }
        
        template <class T, typename rtype> concept xyz_ijk_callable = std::is_invocable<T, ctrs::array<rtype, 3>, int, int, int, int>::value &&
        requires(T t, const ctrs::array<rtype, 3>& x, const int& i, const int& j, const int& k, const int& lb)
        {
            t(x,i,j,k,lb);
        };
        
        template <class T, typename rtype> concept xyz_callable     = std::is_invocable<T, ctrs::array<rtype, 3>>::value                     &&
        requires(T t, const ctrs::array<rtype, 3>& x, const int& i, const int& j, const int& k, const int& lb)
        {
            t(x);
        };
        
        template <class T, typename rtype> concept ijk_callable     = std::is_invocable<T, int, int, int, int>::value                        &&
        requires(T t, const ctrs::array<rtype, 3>& x, const int& i, const int& j, const int& k, const int& lb)
        {
            t(i,j,k,lb);
        };
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <grid::multiblock_grid grid_t, xyz_ijk_callable<typename grid_t::dtype> func_t, ctrs::vec_nd<3, double> vec_t>
        _finline_ auto forward_fillable_args(const grid_t& grid, const func_t& func, const vec_t& x, const int& i, const int& j, const int& k, const int& lb)
        {
            return func(x,i,j,k,lb);
        }
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <grid::multiblock_grid grid_t, xyz_callable<typename grid_t::dtype> func_t, ctrs::vec_nd<3, double> vec_t>
        _finline_ auto forward_fillable_args(const grid_t& grid, const func_t& func, const vec_t& x, const int& i, const int& j, const int& k, const int& lb)
        {
            return func(x);
        }
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <grid::multiblock_grid grid_t, ijk_callable<typename grid_t::dtype> func_t, ctrs::vec_nd<3, double> vec_t>
        _finline_ auto forward_fillable_args(const grid_t& grid, const func_t& func, const vec_t& x, const int& i, const int& j, const int& k, const int& lb)
        {
            return func(i,j,k,lb);
        }
    }
    
    namespace detail
    {
        #pragma GCC diagnostic ignored "-Wattributes"
        template <ctrs::basic_array elem_t, grid::multiblock_array array_t>
        static _finline_ auto unwrap_to_minor_element_type(elem_t& elem, const array_t& arr, const int& i, const int& j, const int& k, const int& lb, const int& maj)
        {
            for (std::size_t n = 0; n < elem.size(); ++n) elem[n] = arr.unwrap_idx(n, i, j, k, lb, maj);
        }
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <typename elem_t, grid::multiblock_array array_t>
        requires (!ctrs::basic_array<elem_t>)
        static _finline_  auto unwrap_to_minor_element_type(elem_t& elem, const array_t& arr, const int& i, const int& j, const int& k, const int& lb, const int& maj)
        {
            elem = arr.unwrap_idx(0, i, j, k, lb, maj);
        }
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <ctrs::basic_array elem_t, grid::multiblock_array array_t>
        static _finline_ auto set_from_minor_element_type(const elem_t& elem, array_t& arr, const int& i, const int& j, const int& k, const int& lb, const int& maj)
        {
            for (std::size_t n = 0; n < elem.size(); ++n)
            {
                arr.unwrap_idx(n, i, j, k, lb, maj) = elem[n];
            }
        }
        
        #pragma GCC diagnostic ignored "-Wattributes"
        template <typename elem_t, grid::multiblock_array array_t>
        requires (!ctrs::basic_array<elem_t>)
        static _finline_  auto set_from_minor_element_type(const elem_t& elem, array_t& arr, const int& i, const int& j, const int& k, const int& lb, const int& maj)
        {
            arr.unwrap_idx(0, i, j, k, lb, maj) = elem;
        }
        
        template <class T> concept has_arg_type = requires(T t) {typename T::arg_type;};
        
        template <typename array_t, typename callable_t> struct converted_elem
        {
            typedef typename array_t::unwrapped_minor_type type;
        };
        
        template <typename array_t, has_arg_type callable_t> struct converted_elem<array_t, callable_t>
        {
            typedef typename callable_t::arg_type type;
        };
        
        //SFINAE for now
        
        //coordinates-only
        template <typename kernel_t, typename array_t, typename idx_t>
        requires (std::invocable<kernel_t, const typename array_t::coord_point_type&>)
        auto invoke_kernel(const kernel_t& kernel, const array_t& arr, const idx_t& idx)
        {
            return kernel(arr.get_grid().get_coords(idx));
        }
        
        template <typename kernel_t, typename array_t, typename idx_t>
        requires (std::invocable<kernel_t, const typename array_t::coord_point_type&, const typename array_t::index_type&>)
        auto invoke_kernel(const kernel_t& kernel, const array_t& arr, const idx_t& idx)
        {
            return kernel(arr.get_grid().get_coords(idx), idx);
        }
        
        template <typename kernel_t, typename array_t, typename idx_t>
        requires (std::invocable<kernel_t, const typename array_t::index_type&>)
        auto invoke_kernel(const kernel_t& kernel, const array_t& arr, const idx_t& idx)
        {
            return kernel(idx);
        }
    }
    
    template <typename array_t, typename lb_idx_t, typename kernel_t>
    void block_elem_loop(array_t& array, const lb_idx_t& lb, const grid::exchange_inclusion_e& exchange_policy, const kernel_t& kernel)
    {
        auto bound_box = array.get_grid_index_bounding_box(exchange_policy);
        const int lb_dim = grid::get_subidx<array.centering_type(), grid::lb_subindex>::value;
        
        //modify the bounding box to only loop over this one block
        bound_box.min(lb_dim) = (int)(lb);
        bound_box.max(lb_dim) = (int)(lb+1);
        
        md_loop(bound_box, kernel);
    }
    
    template <typename data_t, typename array_t, typename elem_t>
    requires(array_t::variable_map_type::rank()==0)
    void set_elem_value(const data_t& data, array_t& arr, const elem_t& elem)
    {
        arr(elem) = data;
    }
    
    template <typename data_t, typename array_t, typename elem_t>
    requires(array_t::variable_map_type::rank()==1)
    void set_elem_value(const data_t& data, array_t& arr, const elem_t& elem)
    {
        for (int i = 0; i < arr.var_map().size() ; ++i) arr(i, elem) = data[i];
    }
    
    template <class array_t, class kernel_t>
    void fill_array(array_t& arr, const kernel_t& kernel, const grid::exchange_inclusion_e& exchange_policy=grid::include_exchanges)
    {
        const auto& grid = arr.get_grid();
        const auto nlb = grid.get_num_local_blocks();
        
        //consider fundamentally separating the block dimension with the ijk dimensions!
        for (auto lb: range(0, nlb))
        {
            const auto loop_func = [&](const auto& elem) -> void
            {
                const auto data = detail::invoke_kernel(kernel, arr, elem);
                set_elem_value(data, arr, elem);
            };
            block_elem_loop(arr, (typename array_t::index_integral_type)lb, exchange_policy, loop_func);
        }    
    }
    
    template <grid::multiblock_array array_t, class callable_t, reduce_ops::reduce_operation<typename array_t::value_type> reduce_t>
    requires std::invocable<callable_t, typename detail::converted_elem<array_t, callable_t>::type>
    auto transform_reduce(
        const array_t& arr,
        const callable_t& func,
        reduce_t& reduce_oper,
        const grid::exchange_inclusion_e& exchange_policy=grid::exclude_exchanges)
    {
        const auto& grid = arr.get_grid();
        auto grid_range = grid.get_range(arr.centering_type(), exchange_policy);
        typename detail::converted_elem<array_t, callable_t>::type init_data;
        detail::unwrap_to_minor_element_type(init_data, arr, 0, 0, 0, 0, 0);
        reduce_oper.init(func(init_data));
        for (auto maj: range(0, arr.get_major_dims().total_size()))
        {
            for (auto i: grid_range)
            {
                typename detail::converted_elem<array_t, callable_t>::type data;
                detail::unwrap_to_minor_element_type(data, arr, i[0], i[1], i[2], i[3], maj);
                reduce_oper.reduce_elem(func(data));
            }
        }
        return grid.group().reduce(reduce_oper.value,reduce_oper.equiv_par_op());
    }
    
    template <grid::multiblock_array array_t, class callable_t>
    requires std::invocable<callable_t, typename detail::converted_elem<array_t, callable_t>::type>
    auto& transform_inplace(array_t& arr, const callable_t& func, const grid::exchange_inclusion_e& exchange_policy=grid::exclude_exchanges)
    {
        const auto& grid = arr.get_grid();
        auto grid_range = grid.get_range(arr.centering_type(), exchange_policy);
        for (auto maj: range(0, arr.get_major_dims().total_size()))
        {
            for (auto i: grid_range)
            {
                typename detail::converted_elem<array_t, callable_t>::type data;
                detail::unwrap_to_minor_element_type(data, arr, i[0], i[1], i[2], i[3], maj);
                detail::set_from_minor_element_type(func(data), arr, i[0], i[1], i[2], i[3], maj);
            }
        }
        return arr;
    }
}
