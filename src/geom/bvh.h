#pragma once

#include <string>
#include <vector>
#include <array>
#include <fstream>

#include "core/bounding_box.h"
#include "core/ctrs.h"

namespace spade::geom
{
    struct bvh_params_t
    {
        int max_level        = 3; //supersedes target_max_elems
        int target_max_elems = 10;
    };
    
    template <const int dim, typename float_t, template <typename> typename container_tt>
    struct bvh_impl_t
    {
        using idx_t    = std::size_t;
        using md_idx_t = ctrs::array<int, dim>;
        static constexpr idx_t no_val = std::string::npos;
        idx_t num_base_squares;
        int degree  = 4;
        using bnd_t  = bound_box_t<float_t,dim>;
        using ibnd_t = bound_box_t<int, dim>;
        using pnt_t  = ctrs::array<float_t, dim>;
        
        constexpr static int bvh_dim() { return dim; }
        
        //idea: we will organize the stored handles by a depth-first traversal
        
        bnd_t glob_bounds;
        container_tt<bnd_t> bounds;
        container_tt<idx_t> data;
        container_tt<idx_t> data_begin;
        container_tt<idx_t> data_end;
        container_tt<idx_t> children;
        container_tt<int>   levels;
        
        bvh_impl_t(){}
        
        void clear_all()
        {
            bounds.clear();
            data.clear();
            data_begin.clear();
            data_end.clear();
            children.clear();
            levels.clear();
        }
        
        idx_t get_max_elems() const
        {
            idx_t output = 0;
            for (idx_t i = 0; i < bounds.size(); ++i) output = utils::max(output, data_end[i] - data_begin[i]);
            return output;   
        }
        
        bnd_t get_subbox(const md_idx_t& ii, const bnd_t& bnd) const
        {
            bnd_t output;
            for (int d = 0; d < bvh_dim(); ++d)
            {
                float_t dx = bnd.size(d)/degree;
                output.min(d) = bnd.min(d) + ii[d]*dx;
                output.max(d) = output.min(d) + dx;
            }
            return output;
        }
        
        template <typename check_t>
        void calculate(
            const bnd_t& bnd,
            const idx_t& lsize,
            const check_t& e_check,
            const bvh_params_t& params = bvh_params_t())
        {
            glob_bounds = bnd;
            this->clear_all();
            num_base_squares = 1;
            for (int i = 0; i < bvh_dim(); ++i) num_base_squares *= degree;
            
            //compute first-level bounding boxes
            for (idx_t i = 0; i < num_base_squares; ++i)
            {
                auto ii = get_index(i);
                bnd_t bnd = get_subbox(ii, glob_bounds);
                bounds.push_back(bnd);
                children.push_back(no_val);
            }
            
            //Now we put every element into the initial bounding boxes
            for (idx_t i = 0; i < num_base_squares; ++i)
            {
                const auto& bnd = bounds[i];
                //This is O(n^2) but in practice it shouldn't be too expensive
                data_begin.push_back(data.size());
                for (idx_t j = 0; j < lsize; ++j)
                {
                    if (e_check(j, bnd))
                    {
                        data.push_back(j);
                    }
                }
                data_end.push_back(data.size());
            }
            
            levels.resize(bounds.size(), 0);
            
            //refine over and over until the limits in params are satisfied
            // this->refine(e_check, params);
            while (this->refine(e_check, params)) {}
        }
        
        template <typename check_t>
        bool refine(const check_t& e_check, const bvh_params_t& params)
        {
            bool did_refine = false;
            const idx_t num_blocks = bounds.size();
            for (idx_t i = 0; i < num_blocks; ++i)
            {
                bool is_terminal       = (children[i] == no_val);
                bool has_too_many_pts  = ((data_end[i] - data_begin[i]) > params.target_max_elems);
                bool less_than_max_lev = (levels[i] < params.max_level);
                if (is_terminal && has_too_many_pts && less_than_max_lev)
                {
                    did_refine  = true;
                    children[i] = bounds.size();
                    
                    //Refining the current node
                    for (idx_t j = 0; j < num_base_squares; ++j)
                    {
                        bnd_t bnd = get_subbox(get_index(j), bounds[i]);
                        bounds.push_back(bnd);
                        children.push_back(no_val);
                        levels.push_back(levels[i]+1);
                        
                        data_begin.push_back(data.size());
                        for (idx_t p = data_begin[i]; p < data_end[i]; ++p)
                        {
                            if (e_check(data[p], bnd))
                            {
                                data.push_back(data[p]);
                            }
                        }
                        data_end.push_back(data.size());
                    }
                }
            }
            return did_refine;
        }
        
        ibnd_t get_index_range(const bnd_t& bvh_bnd, const bnd_t& bnd) const
        {
            ibnd_t output;
            for (int d = 0; d < bvh_dim(); ++d)
            {
                const auto dx = bvh_bnd.size(d) / degree;
                
                output.min(d) = ((bnd.min(d)-bvh_bnd.min(d))/dx);
                output.max(d) = ((bnd.max(d)-bvh_bnd.min(d))/dx) + 1;
                
                output.min(d) = utils::max(output.min(d), 0);
                output.max(d) = utils::max(output.max(d), 0);
                
                output.min(d) = utils::min(output.min(d), degree);
                output.max(d) = utils::min(output.max(d), degree);
            }
            return output;
        }
        
        md_idx_t get_index_range(const bnd_t& bvh_bnd, const pnt_t& pnt) const
        {
            md_idx_t output;
            for (int d = 0; d < bvh_dim(); ++d)
            {
                const auto dx = bvh_bnd.size(d) / degree;
                output[d]     = (pnt[d]-bvh_bnd.min(d))/dx;
            }
            return output;
        }
        
        template <typename iter_t>
        void check_elements(const iter_t& iter, const bnd_t& bnd, const bnd_t& bvh_bnd, const idx_t& patch_idx) const
        {
            auto ibnd = get_index_range(bvh_bnd, bnd);
            const auto loop = [&](const auto& ii)
            {
                idx_t idx = patch_idx+get_index(ii);
                // print(idx, children.size());
                if (children[idx] == no_val)
                {
                    //Terminal node: loop through all the elements in this region
                    for (idx_t i = data_begin[idx]; i < data_end[idx]; ++i)
                    {
                        iter(data[i]);
                    }
                }
                else
                {
                    const auto& next_bvh_bnd = bounds[idx];
                    check_elements(iter, bnd, next_bvh_bnd, children[idx]);
                }
            };
            
            dispatch::execute(ibnd, loop, device::cpu);
        }
        
        template <typename iter_t>
        void check_elements(const iter_t& iter, const pnt_t& pnt, const bnd_t& bvh_bnd, const idx_t& patch_idx) const
        {
            if (!bvh_bnd.contains(pnt)) return;
            const auto md      = get_index_range(bvh_bnd, pnt);
            const auto idx     = patch_idx + get_index(md);
            if (children[idx] == no_val)
            {
                for (idx_t i = data_begin[idx]; i < data_end[idx]; ++i)
                {
                    iter(data[i]);
                }
            }
            else
            {
                check_elements(iter, pnt, bounds[idx], children[idx]);
            }
        }
        
        template <typename iter_t>
        void check_elements(const iter_t& iter, const bnd_t& bnd) const
        {
            check_elements(iter, bnd, glob_bounds, 0);
        }
        
        template <typename iter_t>
        void check_elements(const iter_t& iter, const pnt_t& pnt) const
        {
            check_elements(iter, pnt, glob_bounds, 0);
        }
        
        // Rare pass by value!!
        md_idx_t get_index(idx_t i) const
        {
            md_idx_t output;
            for (int j = 0; j < output.size(); ++j)
            {
                output[j] = i % degree;
                i -= output[j];
                i /= degree;
            }
            return output;
        }
        
        idx_t get_index(const md_idx_t& i) const
        {
            idx_t output = 0;
            idx_t coeff = 1;
            for (int j = 0; j < i.size(); ++j)
            {
                output += coeff*i[j];
                coeff  *= degree;
            }
            return output;
        }
    };
    
    template <const int dim, typename float_t> using bvh_t = bvh_impl_t<dim, float_t, std::vector>;
    
    namespace detail
    {
        static void debug_output_bvh(const std::string& fname, const auto& bvh)
        {
            static_assert((bvh.bvh_dim() == 2) || (bvh.bvh_dim() == 3), "can only debug 2d or 3d bvh");
            std::size_t npts   = 0;
            std::size_t nboxes = 0;
            for (std::size_t i = 0; i < bvh.bounds.size(); ++i)
            {
                if (bvh.children[i] == bvh.no_val)
                {
                    nboxes++;
                    npts += 8;
                }
            }
            
            std::ofstream mf(fname);
            mf << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << npts << " double\n";
            for (std::size_t i = 0; i < bvh.bounds.size(); ++i)
            {
                if (bvh.children[i] == bvh.no_val)
                {
                    const auto& bnd = bvh.bounds[i];
                    if constexpr (bvh.bvh_dim() == 3)
                    {
                        mf << bnd.min(0) << " " << bnd.min(1) << " " << bnd.min(2) << "\n";
                        mf << bnd.max(0) << " " << bnd.min(1) << " " << bnd.min(2) << "\n";
                        mf << bnd.min(0) << " " << bnd.max(1) << " " << bnd.min(2) << "\n";
                        mf << bnd.max(0) << " " << bnd.max(1) << " " << bnd.min(2) << "\n";
                        mf << bnd.min(0) << " " << bnd.min(1) << " " << bnd.max(2) << "\n";
                        mf << bnd.max(0) << " " << bnd.min(1) << " " << bnd.max(2) << "\n";
                        mf << bnd.min(0) << " " << bnd.max(1) << " " << bnd.max(2) << "\n";
                        mf << bnd.max(0) << " " << bnd.max(1) << " " << bnd.max(2) << "\n";
                    }
                    else
                    {
                        mf << bnd.min(0) << " " << bnd.min(1) << " " << 0.0 << "\n";
                        mf << bnd.max(0) << " " << bnd.min(1) << " " << 0.0 << "\n";
                        mf << bnd.min(0) << " " << bnd.max(1) << " " << 0.0 << "\n";
                        mf << bnd.max(0) << " " << bnd.max(1) << " " << 0.0 << "\n";
                        mf << bnd.min(0) << " " << bnd.min(1) << " " << 1.0 << "\n";
                        mf << bnd.max(0) << " " << bnd.min(1) << " " << 1.0 << "\n";
                        mf << bnd.min(0) << " " << bnd.max(1) << " " << 1.0 << "\n";
                        mf << bnd.max(0) << " " << bnd.max(1) << " " << 1.0 << "\n";
                    }
                }
            }
            
            mf << "CELLS " << nboxes << " " << 9*nboxes << "\n";
            
            std::size_t ipt = 0;
            for (std::size_t i = 0; i < bvh.bounds.size(); ++i)
            {
                if (bvh.children[i] == bvh.no_val)
                {
                    mf << "8";
                    for (int pp = 0; pp < 8; ++pp)
                    {
                        mf << " " << ipt++;
                    }
                    mf << "\n";
                }
            }
            
            mf << "CELL_TYPES " << nboxes << "\n";
            for (std::size_t pp = 0; pp < nboxes; ++pp)
            {
                mf << 11 << "\n";
            }
            
            mf << "CELL_DATA " << nboxes << "\n";
            
            const auto write_scalar = [&](const std::string& scname, const auto& get_data)
            {
                mf << "SCALARS " << scname << " double\n";
                mf << "LOOKUP_TABLE default\n";
                for (std::size_t i = 0; i < bvh.bounds.size(); ++i)
                {
                    if (bvh.children[i] == bvh.no_val)
                    {
                        mf << get_data(i) << "\n";
                    }
                }
            };
            
            write_scalar("count", [&](const auto i) { return bvh.data_end[i] - bvh.data_begin[i]; });
        }
    }
}