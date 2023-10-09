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
        int max_level        = 6; //supersedes max_elems
        int target_max_elems = 10;
    };
    
    template <const int dim, typename float_t, template <typename> typename container_tt = std::vector>
    struct bvh_t
    {
        using idx_t    = std::size_t;
        using md_idx_t = ctrs::array<idx_t, dim>;
        static constexpr idx_t no_val = std::string::npos;
        int degree  = 4;
        using bnd_t = bound_box_t<float_t,dim>;
        
        constexpr static int bvh_dim() { return dim; }
        
        //idea: we will organize the stored handles by a depth-first traversal
        
        bnd_t glob_bounds;
        container_tt<bnd_t> bounds;
        container_tt<idx_t> data;
        container_tt<idx_t> data_begin;
        container_tt<idx_t> data_end;
        container_tt<idx_t> children;
        container_tt<idx_t> parents;
        
        bvh_t(){}
        
        template <typename element_check_t, typename element_bbox_t>
        void calculate(
            const bnd_t& bnd,
            const idx_t& lsize,
            const element_check_t& e_check,
            const element_bbox_t&  e_bbx,
            const bvh_params_t& params = bvh_params_t())
        {
            glob_bounds = bnd;
            data.clear();
            children.clear();
            parents.clear();
            idx_t num_base_squares = 1;
            for (int i = 0; i < bvh_dim(); ++i) num_base_squares *= degree;
            
            //compute first-level bounding boxes
            for (idx_t i = 0; i < num_base_squares; ++i)
            {
                auto ii = get_index(i);
                bnd_t bnd;
                for (int d = 0; d < bvh_dim(); ++d)
                {
                    float_t dx = glob_bounds.size(d)/degree;
                    bnd.min(d) = glob_bounds.min(d) + ii[d]*dx;
                    bnd.max(d) = bnd.min(d) + dx;
                }
                bounds.push_back(bnd);
                children.push_back(no_val);
                parents.push_back(no_val);
            }
            
            //Now we put every element into the initial bounding boxes
            for (idx_t i = 0; i < num_base_squares; ++i)
            {
                const auto& bnd = bounds[i];
                //This is O(n^2) but in practice it shouldn't be too expensive
                data_begin.push_back(data.size());
                for (idx_t j = 0; j < lsize; ++j)
                {
                    if (e_check(i, bnd))
                    {
                        data.push_back(i);
                    }
                }
                data_end.push_back(data.size());
            }
        }
        
        // idx = i0 + n*i1 + n*n*i2 + n*n*n*i3
        
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
        }
    }
}