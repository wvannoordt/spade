#pragma once

#include <tuple>
#include "grid/grid.h"
#include "algs/transform_inplace.h"

namespace spade::sampling
{
    template <typename grid_t>
    struct slice_params_t
    {
        using grid_type   = grid_t;
        using slice_type  = typename grid_type::slice_type;
        using coord_type  = typename grid_type::coord_type;
        using tagged_type = decltype(utils::tag[partition::global](0UL));
        coord_type x_slice;
        int dir_slice;
        ctrs::array<int, 2> tan_dirs;
        bool is_collapsible;
        std::vector<tagged_type> source_block;   // block number on 2d grid to source block on 3d grid
        std::vector<tagged_type> collapse_block; // block number on 3d grid to dest.  block on 2d grid, empty if not collapsible
        
        bool collapsible() const
        {
            return is_collapsible;
        }
        
        template <typename iarr_t>
        _sp_hybrid auto reduce_array(const iarr_t& iarr) const
        {
            using output_t = typename iarr_t::template change_size<2>;
            output_t output;
            output[0] = iarr[tan_dirs[0]];
            output[1] = iarr[tan_dirs[1]];
            return output;
        }
    };
    
    template <typename grid_t>
    inline [[nodiscard]] auto slice_grid(const grid_t& grid, const typename grid_t::coord_type& x_slice, int dir)
    {
        using ogrid_t  = typename grid_t::slice_type;
        using coord_t  = typename grid_t::coord_type;
        
        if (dir < 0 || dir >= grid_t::dim())
        {
            std::string msg = "Invalid direction argument to slice_grid: expecting value between 0 and " + std::to_string(grid_t::dim()-1);
            msg += ", but was passed " + std::to_string(dir);
            throw except::sp_exception(msg);
        }
        
        slice_params_t<grid_t> params_output;
        params_output.x_slice   = x_slice;
        params_output.dir_slice = dir;
        int d0 = (dir + 1) % grid_t::dim();
        int d1 = (dir + 2) % grid_t::dim();
        params_output.tan_dirs = ctrs::make_array(utils::min(d0, d1), utils::max(d0, d1));
        
        using oblocks_t  = typename ogrid_t::blocks_type;
        using desig_type = typename oblocks_t::desig_type;
        
        const auto& blocks = grid.get_blocks();
        
        using bbx_t = bound_box_t<coord_t, desig_type::size()>;
        bbx_t obounds;
        for (int i = 0; i < 2; ++i)
        {
            obounds.min(i) = blocks.bounds.min(params_output.tan_dirs[i]);
            obounds.max(i) = blocks.bounds.max(params_output.tan_dirs[i]);
        }
        
        oblocks_t oblocks(params_output.reduce_array(blocks.num_blocks), obounds);
        ogrid_t  ogrid(params_output.reduce_array(grid.get_num_cells()), oblocks, grid.coord_sys(), grid.group());
        
        using lbglob_t = decltype(utils::tag[partition::global](0UL));
        using ref_t    = oblocks_t::refine_type;
        
        const coord_t eps = 0.01*grid.compute_dx_min()[dir];
        auto xx0 = blocks.bounds.min(dir);
        auto xx1 = blocks.bounds.max(dir);
        bool close_to_border = (utils::abs(x_slice - xx0) < eps) || (utils::abs(x_slice - xx1) < eps);
        
        if (x_slice < xx0 || x_slice > xx1)
        {
            std::stringstream oss0, oss1;
            oss0 << x_slice;
            oss1 << blocks.bounds;
            std::string mesg = "Attempted to slice a grid with bad coordinate in direction ";
            mesg += std::to_string(dir);
            mesg += ": coordinate " + oss0.str() + " is outside bounding box " + oss1.str();
            throw except::sp_exception(mesg);
        }
        
        const auto is_sliced = [&](const auto& lbglob)
        {
            const auto bbx = grid.get_bounding_box(lbglob);
            auto x0 = bbx.min(dir);
            auto x1 = bbx.max(dir);
            
            // Awful implementation
            if (close_to_border)
            {
                return (x_slice >= x0) && (x_slice <= x1);
            }
            else
            {
                bool close_to_block_boundary = (utils::abs(x_slice - x0) < eps) || (utils::abs(x_slice - x1) < eps);
                if (close_to_block_boundary)
                {
                    x0 += 0.5*eps;
                    x1 += 0.5*eps;
                    return (x_slice >= x0) && (x_slice <= x1);
                }
                else
                {
                    return (x_slice >= x0) && (x_slice <= x1);
                }
            }
            return false;
        };
        
        auto slice_blocks = grid.select_blocks(is_sliced, partition::global);
        
        geom::bvh_t<ogrid_t::dim(), coord_t> block_bvh;
        geom::bvh_params_t bvhparams{3, 1000};
        
        const auto el_check = [&](const std::size_t& vec_idx, const auto& bnd_in)
        {
            const auto lb = slice_blocks[vec_idx];
            auto block_box = grid.get_bounding_box(lb);
            const auto dx = grid.get_dx(lb);
            bound_box_t<coord_t, ogrid_t::dim()> block_box_2d;
            for (int d = 0; d < ogrid_t::dim(); ++d)
            {
                block_box_2d.min(d) = block_box.min(params_output.tan_dirs[d]) - 2*dx[params_output.tan_dirs[d]];
                block_box_2d.max(d) = block_box.max(params_output.tan_dirs[d]) + 2*dx[params_output.tan_dirs[d]];
            }
            return block_box_2d.intersects(bnd_in);
        };
        
        bound_box_t<coord_t, ogrid_t::dim()> bnd;
        for (int d = 0; d < ogrid_t::dim(); ++d)
        {
            bnd.min(d) = blocks.bounds.min(params_output.tan_dirs[d]);
            bnd.max(d) = blocks.bounds.max(params_output.tan_dirs[d]);
        }
        block_bvh.calculate(bnd, slice_blocks.size(), el_check, bvhparams);
        
        bool refine_done = false;
        while (!refine_done)
        {
            std::vector<lbglob_t> to_refine;
            std::vector<ref_t>    ref_dirs;
            auto dxm = ogrid.compute_dx_min();
            const coord_t dx_search = utils::min(dxm[0], dxm[1]);
            for (std::size_t ilb = 0; ilb < ogrid.get_num_global_blocks(); ++ilb)
            {
                lbglob_t lb_glob  = utils::tag[partition::global](ilb);
                ref_t    all_dirs = false;
                ref_t    all_dirs_lo = false;
                auto ibbx = ogrid.get_bounding_box(lb_glob);
                auto xc   = ibbx.center();
                
                bound_box_t<coord_t, ogrid_t::dim()> ibbx2d;
                ibbx2d.min(0) = ibbx.min(params_output.tan_dirs[0]);
                ibbx2d.max(0) = ibbx.max(params_output.tan_dirs[0]);
                ibbx2d.min(1) = ibbx.min(params_output.tan_dirs[1]);
                ibbx2d.max(1) = ibbx.max(params_output.tan_dirs[1]);
                ctrs::array<coord_t, 2> xc2d{xc[0], xc[1]};
                ctrs::array<coord_t, 2> bbx_dx{dx_search, dx_search};
                auto bnd = utils::bbox_around(xc2d, bbx_dx);
                const auto dx_2d = ogrid.get_dx(lb_glob);
                xc[dir] = x_slice;
                ctrs::array<ctrs::array<bool, 3>,2> special_case_check;
                special_case_check[0] = false;
                special_case_check[1] = false;
                const auto eval = [&](const std::size_t& vec_idx_cand)
                {
                    auto lb_tagged = slice_blocks[vec_idx_cand];
                    auto bbx_igrd  = grid.get_bounding_box(lb_tagged);
                    
                    bound_box_t<coord_t, ogrid_t::dim()> vbbx2d;
                    vbbx2d.min(0) = bbx_igrd.min(params_output.tan_dirs[0]);
                    vbbx2d.max(0) = bbx_igrd.max(params_output.tan_dirs[0]);
                    vbbx2d.min(1) = bbx_igrd.min(params_output.tan_dirs[1]);
                    vbbx2d.max(1) = bbx_igrd.max(params_output.tan_dirs[1]);
                    
                    const auto ixc = ibbx2d.center();
                    const auto oxc = vbbx2d.center();
                    if (ibbx2d.contains(oxc))
                    {
                        
                        const auto dx_3d = grid.get_dx(lb_tagged);
                        for (int d = 0; d < 2; ++d)
                        {
                            const auto dx_d3_i = dx_3d[params_output.tan_dirs[d]];
                            int icase = -1;
                            if (utils::abs(ixc[d] - oxc[d]) < 0.2*dx_d3_i) icase = 1;
                            else if (ixc[d] < oxc[d]) icase = 0;
                            else if (ixc[d] > oxc[d]) icase = 2;
                            special_case_check[d][icase] = true;
                            all_dirs[d] = all_dirs[d] || ((dx_2d[d] > 1.01*dx_d3_i));
                        }
                    }
                };
                block_bvh.check_elements(eval, bnd);
                
                bool exclude0 = special_case_check[0][0] && special_case_check[0][1] && special_case_check[0][2];
                bool exclude1 = special_case_check[1][0] && special_case_check[1][1] && special_case_check[1][2];
                
                all_dirs[0] = all_dirs[0] && !exclude0;
                all_dirs[1] = all_dirs[1] && !exclude1;
                if (all_dirs[0] || all_dirs[1])
                {
                    to_refine.push_back(lb_glob);
                    ref_dirs.push_back(all_dirs);
                }
            }
            
            refine_done = (to_refine.size() == 0);
            
            if (!refine_done)
            {
                ctrs::array<bool, 2> peri = false; //Note: no need for periodic refinement or constraints here!
                ogrid.refine_blocks(to_refine, peri, ref_dirs, amr::constraints::none);
            }
        }
        params_output.source_block.resize(ogrid.get_num_global_blocks(), ogrid.get_num_global_blocks() + 1);
        params_output.collapse_block.resize(grid.get_num_global_blocks(), grid.get_num_global_blocks() + 1);
        
        const auto ogrid_el_check = [&](const std::size_t& lb_ogrid, const auto& bnd_in)
        {
            const auto lb = utils::tag[partition::global](lb_ogrid);
            auto block_box = ogrid.get_bounding_box(lb);
            const auto dx = ogrid.get_dx(lb);
            bound_box_t<coord_t, ogrid_t::dim()> block_box_2d;
            for (int d = 0; d < ogrid_t::dim(); ++d)
            {
                block_box_2d.min(d) = block_box.min(d) - 2*dx[d];
                block_box_2d.max(d) = block_box.max(d) + 2*dx[d];
            }
            block_box_2d = block_box_2d.inflate(coord_t(1.02));
            return block_box_2d.intersects(bnd_in.inflate(1.02));
        };
        block_bvh.calculate(obounds, ogrid.get_num_global_blocks(), ogrid_el_check, bvhparams);
        
        const auto dxmin = ogrid.compute_dx_min();
        const auto eps2 = 0.01*utils::min(dxmin[0], dxmin[1]);
        for (std::size_t ilb = 0; ilb < grid.get_num_global_blocks(); ++ilb)
        {
            const auto ilb_tagged = utils::tag[partition::global](ilb);
            bound_box_t<coord_t, ogrid_t::dim()> grid_bbx_2d;
            
            const auto bbx = grid.get_bounding_box(ilb_tagged);
            grid_bbx_2d.min(0) = bbx.min(params_output.tan_dirs[0]);
            grid_bbx_2d.max(0) = bbx.max(params_output.tan_dirs[0]);
            grid_bbx_2d.min(1) = bbx.min(params_output.tan_dirs[1]);
            grid_bbx_2d.max(1) = bbx.max(params_output.tan_dirs[1]);
            
            const auto xc_3d_grid = grid_bbx_2d.center();
            
            const auto check = [&](const std::size_t& lb_ogrid)
            {
                auto lb_tagged = utils::tag[partition::global](lb_ogrid);
                auto bbx_igrd  = ogrid.get_bounding_box(lb_tagged);
                
                bound_box_t<coord_t, ogrid_t::dim()> vbbx2d;
                vbbx2d.min(0) = bbx_igrd.min(0);
                vbbx2d.max(0) = bbx_igrd.max(0);
                vbbx2d.min(1) = bbx_igrd.min(1);
                vbbx2d.max(1) = bbx_igrd.max(1);
                const auto test_xc_2d_grid = vbbx2d.center();
                const auto diff = ctrs::array_norm(xc_3d_grid - test_xc_2d_grid);
                if (diff < eps2)
                {
                    params_output.collapse_block[ilb_tagged.value] = lb_tagged;
                    
                    if (is_sliced(ilb_tagged))
                    {
                        params_output.source_block[lb_tagged.value] = ilb_tagged;
                    }
                }
            };
            
            block_bvh.check_elements(check, grid_bbx_2d);
        }
        
        for (const auto& lb: params_output.source_block)
        {
            if (lb.value >= grid.get_num_global_blocks())
            {
                throw except::sp_exception("Error after performing grid slice: cannot find source block for at least one slice block.");
            }
        }
        
        params_output.is_collapsible = true;
        for (const auto& lb: params_output.collapse_block)
        {
            if (lb.value >= ogrid.get_num_global_blocks())
            {
                params_output.is_collapsible = false;
            }
        }
        
        if (!params_output.is_collapsible)
        {
            params_output.collapse_block.clear();
        }
        
        return std::make_tuple(ogrid, params_output);
    }
    
    namespace detail
    {
        template <typename sarr_t, typename arr_t, typename params_t, typename kern_t>
        inline void slice_arr_impl(sarr_t& sarr, const arr_t& arr, const params_t& params, const kern_t& kern)
        {
            using coord_t = typename arr_t::grid_type::coord_type;
            const auto& sarr_grid = sarr.get_grid();
            const auto& arr_grid  = arr.get_grid();            
            
            if (sarr_grid.group().size() > 1)
            {
                throw except::sp_exception("slice_array is currently unsupported for multi-gpu implementation");
            }
            
            const auto sarr_grid_img = sarr_grid.image(partition::local, sarr.device());
            const auto arr_grid_img  = arr_grid.image(partition::local, arr.device());
            
            auto sarr_img = sarr.image();
            auto arr_img  = arr.image();
            
            device::shared_vector<std::size_t> lookup;
            lookup.resize(params.source_block.size());
            for (std::size_t ilb = 0; ilb < params.source_block.size(); ++ilb)
            {
                lookup[ilb] = params.source_block[ilb].value;
            }
            lookup.transfer();
            
            const auto lkp_img = utils::make_vec_image(lookup.data(sarr.device()));
            
            int idir           = params.dir_slice;
            const auto x_slice = params.x_slice;
            int d0 = params.tan_dirs[0];
            int d1 = params.tan_dirs[1];
            const auto tform = [=] _sp_hybrid (const typename sarr_t::index_type& ii)
            {
                auto iii        = ii;
                iii.lb()        = lkp_img[ii.lb()];
                iii.i(d0)       = ii.i();
                iii.i(d1)       = ii.j();
                const auto dx_i = arr_grid_img.get_dx(idir, iii.lb());
                const auto bbx  = arr_grid_img.get_bounding_box(iii.lb());                
                int i_norm      = int(((x_slice - bbx.min(idir)) / dx_i) - coord_t(0.5));
                iii.i(idir)     = i_norm;
                const auto data = algs::invoke_at(arr_grid_img, arr_img, iii, kern);
                return data;
            };
            
            algs::transform_inplace(sarr, tform);
        }
    }
    
    template <typename array_t, typename sgrid_t, typename func_t>
    static auto slice_array(const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params, const func_t& func)
    {
        const auto kern = omni::to_omni<array_t::centering_type()>(func, array);
        using kern_t = typename utils::remove_all<decltype(kern)>::type;
        using new_alias_type = typename kern_t::output_type;
        using output_type = typename array_t::slice_type<new_alias_type>;
        
        new_alias_type fl;
        const auto nexch = params.reduce_array(array.get_num_exchange());
        // output_type output(sgrid, fl, nexch, array.device(), array.mmap_tag);
        output_type output(sgrid, fl, nexch, array.device(), mem_map::linear);
        
        output = 0;
        
        detail::slice_arr_impl(output, array, params, func);
        
        return output;
    }
    
    template <typename array_t, typename sgrid_t>
    inline typename array_t::slice_type<typename array_t::alias_type>
    slice_array(
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params)
    {
        const auto lam = [] _sp_hybrid (const typename array_t::alias_type& q) { return q; };
        const auto kern = omni::to_omni<array_t::centering_type()>(lam, array);
        return slice_array(array, sgrid, params, kern);
    }
    
    template <typename alias_t, typename array_t, typename sgrid_t, typename func_t>
    inline void slice_array(
        typename array_t::slice_type<alias_t>& sarray,
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params, const func_t& func)
    {
        detail::slice_arr_impl(sarray, array, params, func);
    }
    
    template <typename array_t, typename sgrid_t>
    inline void slice_array(
        typename array_t::slice_type<typename array_t::alias_type>& sarray,
        const array_t& array, const sgrid_t& sgrid, const slice_params_t<typename array_t::grid_type>& params)
    {
        const auto lam = [] _sp_hybrid (const typename array_t::alias_type& q) { return q; };
        const auto kern = omni::to_omni<array_t::centering_type()>(lam, array);
        detail::slice_arr_impl(sarray, array, params, kern);
    }
}