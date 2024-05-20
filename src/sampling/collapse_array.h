#pragma once

#include <concepts>
#include "sampling/array_slice.h"

namespace spade::sampling
{
    struct id_pair_t { std::size_t lb2d, lb3d; };
    template <typename sarray_t, typename array_t, typename sgrid_t, typename func_t>
    static void collapse_array(
        sarray_t& sarray,
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params, const func_t& func)
    {
        const auto kern = omni::to_omni<array_t::centering_type()>(func, array);
        
        using fund_t = typename array_t::value_type;
        static_assert(std::floating_point<fund_t>, "cannot collapse array that doesn not contain floating-point values");
        if (array.get_grid().group().size() > 1) throw except::sp_exception("collapse_array not yet supported for multi-threading");
        if (!params.collapsible())               throw except::sp_exception("attempted to collapse array with un-collapsible grid");
        
        device::shared_vector<id_pair_t> table;
        table.resize(params.collapse_block.size());
        for (std::size_t i = 0; i < params.collapse_block.size(); ++i)
        {
            table[i] = id_pair_t{params.collapse_block[i], i};
        }
        std::sort(table.begin(), table.end(), [](const auto& a, const auto& b){ return a.lb2d < b.lb2d; });
        table.transfer();
        device::shared_vector<std::size_t> denominator2d;
        denominator2d.resize(sgrid.get_num_global_blocks());
        for (std::size_t i = 0; i < table.size(); ++i)
        {
            denominator2d[table[i].lb2d]++;
        }
        denominator2d.transfer();
        std::size_t max_denom = 0;
        for (const auto& d: denominator2d) max_denom = utils::max(max_denom, d);
        const std::size_t invalid = array.get_grid().get_num_global_blocks() + 1;
        device::shared_vector<std::size_t> list_start;
        device::shared_vector<std::size_t> list_size;
        
        list_start.resize(sgrid.get_num_global_blocks(), invalid);
        list_size.resize(sgrid.get_num_global_blocks(), 0);
        
        for (std::size_t i = 0; i < table.size(); ++i)
        {
            std::size_t lb2d = table[i].lb2d;
            if (list_start[lb2d] == invalid) list_start[lb2d] = i;
            list_size[lb2d]++;
        }
        list_start.transfer();
        list_size.transfer();
        
        using index_type = typename array_t::index_type;
        using val_type   = typename sarray_t::alias_type;
        
        const auto arr_device     = array.device();
        const auto grid_img       = array.get_grid().image(partition::global, arr_device);
        const auto arr_img        = array.image();
        const auto list_start_img = utils::make_vec_image(list_start.data(arr_device));
        const auto list_size_img  = utils::make_vec_image(list_size.data(arr_device));
        const auto table_img      = utils::make_vec_image(table.data(arr_device));
        const int idir   = params.dir_slice;
        const int nx_dir = array.get_grid().get_num_cells(idir);
        const int d0     = params.tan_dirs[0];
        const int d1     = params.tan_dirs[1];
        sarray = fund_t(0.0);
        for (int i_block = 0; i_block < max_denom; ++i_block)
        {
            const auto lam = [=] _sp_hybrid (const index_type& i2d, const val_type& q2d)
            {
                val_type q_out = q2d;
                int lb2d  = i2d.lb();
                int idx   = list_start_img[lb2d];
                int lsize = list_size_img[lb2d];
                idx += (i_block % lsize);
                int lb3d = table_img[idx].lb3d;
                
                index_type i3d;
                i3d.lb()    = lb3d;
                i3d.i(d0)   = i2d.i();
                i3d.i(d1)   = i2d.j();
                i3d.i(idir) = 0;
                
                for (int ii = 0; ii < nx_dir; ++ii)
                {
                    i3d.i(idir) = ii;
                    const auto vl = algs::invoke_at(grid_img, arr_img, i3d, kern);
                    q_out += vl;
                }
                return q_out;
            };
            algs::fill_array(sarray, lam, grid::exclude_exchanges);
        }
        const fund_t coeff = fund_t(1.0)/(max_denom*nx_dir);
        sarray *= coeff;
    }
    
    template <typename array_t, typename sgrid_t>
    static typename array_t::slice_type<typename array_t::alias_type>
    collapse_array(
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params)
    {
        using output_t  = typename array_t::slice_type<typename array_t::alias_type>;
        using alias_t   = typename array_t::alias_type;
        const auto lam  = [=] _sp_hybrid (const alias_t& q) { return q; };
        output_t sarray = slice_array(array, sgrid, params, lam);
        collapse_array(sarray, array, sgrid, params, lam);
        return sarray;
    }
    
    template <typename array_t, typename sgrid_t, typename lam_t>
    static auto
    collapse_array(
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params, const lam_t& lam)
    {
        auto sarray = slice_array(array, sgrid, params, lam);
        collapse_array(sarray, array, sgrid, params, lam);
        return sarray;
    }
    
    template <typename array_t, typename sgrid_t>
    static void collapse_array(
        typename array_t::slice_type<typename array_t::alias_type>& sarray,
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params)
    {
        collapse_array(sarray, array, sgrid, params, [=] _sp_hybrid (const typename array_t::alias_type& q) { return q; });
    }
}