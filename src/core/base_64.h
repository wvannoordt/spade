#pragma once

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <errno.h>
#include <cstring>
#include "core/print.h"
#include "grid/grid.h"
#include "dispatch/support_of.h"
#include "dispatch/execute.h"

namespace spade::io::detail
{
    
    template <typename data_t, typename idx_t>
    requires (ctrs::basic_array<data_t>)
    _sp_hybrid static auto ith_elem(const data_t& d, const idx_t& i) { return d[i]; }
    
    template <typename data_t, typename idx_t>
    _sp_hybrid static auto ith_elem(const data_t& d, const idx_t&) { return d; }
        
    template <class data_t> static inline void stream_base_64(std::ostream& strm, const data_t* data, const size_t& data_size)
    {
        //This is a very slow implentation for now.
        const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        char term = '=';
        char b0, b1, b2, c0, c1, c2, c3;
        std::size_t size = data_size*sizeof(data_t);
        size_t numWriteTotal = (size - size%3)/3;
        char* bytes = (char*)data;
        for (size_t pos = 0; pos < numWriteTotal; pos++)
        {
            b0 = bytes[3*pos];
            b1 = bytes[3*pos+1];
            b2 = bytes[3*pos+2];
            c0 = ((b0 & 0b11111100) >> 2);
            c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
            c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
            c3 = ((b2 & 0b00111111));
            strm << table[c0] << table[c1] << table[c2] << table[c3];
        }
        int numLeft = size%3;
        if (numLeft==0) return;
        char last3[3] = {0};
        char last4[4] = {0};
        int num2write[3] = {0, 2, 3};
        int numTerms[3] = {0, 2, 1};
        for (int i = 0; i < numLeft; i++)
        {
            last3[i] = bytes[size-(size%3)+i];
        }
        b0 = last3[0];
        b1 = last3[1];
        b2 = last3[2];
        c0 = ((b0 & 0b11111100) >> 2);
        c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
        c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
        c3 = ((b2 & 0b00111111));
        last4[0] = c0;
        last4[1] = c1;
        last4[2] = c2;
        last4[3] = c3;
        for (int i = 0; i < num2write[numLeft]; i++)
        {
            strm << table[last4[i]];
        }
        for (int i = 0; i < numTerms[numLeft]; i++)
        {
            strm << term;
        }
    }
    
    template<typename data_t, typename alloc_t>
    static inline void stream_base_64(std::ostream& strm, const std::vector<data_t, alloc_t>& data)
    {
        unsigned int bytes = data.size() * sizeof(data_t);
        stream_base_64(strm, &bytes, 1);
        stream_base_64(strm, data.data(), data.size());
    }
    
    template <typename vchars_t, typename vdata_t, typename vcoord_t, typename grid_t, typename arr_t, typename lbglob_t>
    static inline std::pair<std::size_t, std::size_t>
    vtk_b64_dev_impl(vchars_t& obuf, vdata_t& data_raw, vcoord_t& coord_raw, const grid_t& grid, const arr_t& arr, const lbglob_t& lb_glob)
    {
        timing::tmr_t otmr;
        // timing::scoped_tmr_t stmr("blk");
        using alias_type = arr_t::alias_type;
        constexpr int num_vars = [&]()
        {
            if constexpr (ctrs::basic_array<alias_type>) return alias_type::size();
            else return 1;
        }();
        auto& target = obuf.data(arr.device());
        const auto g_img = grid.image(partition::global, arr.device());
        auto range = dispatch::support_of(arr, grid::exclude_exchanges);
        auto range2 = range;
        range2.upper.i()++;
        range2.upper.j()++;
        range2.upper.k()++;
        int the_lb      = grid.get_partition().to_local(lb_glob).value;
        int the_lb_glob = lb_glob.value;
        range.lower.lb() = 0;
        range.upper.lb() = num_vars;
        
        range2.lower.lb() = 0;
        range2.upper.lb() = 1;
        
        std::size_t ni = grid.get_num_cells(0);
        std::size_t nj = grid.get_num_cells(1);
        std::size_t nk = grid.get_num_cells(2);
        
        std::size_t nni = ni+1;
        std::size_t nnj = nj+1;
        std::size_t nnk = nk+1;
        
        auto nx = ctrs::make_array(ni, nj, nk);
        
        std::size_t voffst = ni*nj*nk;
        
        using data_t     = vdata_t::value_type;
        using coor_t     = vcoord_t::value_type;
        
        auto& vdr = data_raw.data(arr.device());
        auto& cdr = coord_raw.data(arr.device());
        
        data_t* ptr  = &vdr[0];
        coor_t* cptr = &cdr[0];
        auto img = arr.image();
        using index_type = decltype(range)::index_type;
        auto load0 = [=] _sp_hybrid (const index_type& iii) mutable
        {
            auto i = iii;
            int iv = i.lb();
            i.lb() = the_lb;
            
            int v = iv;
            std::size_t offst = i.i() + ni*i.j() + ni*nj*i.k() + voffst*v;
            const auto data   = img.get_elem(i);
            ptr[offst] = ith_elem(data, v);
        };
        
        auto load1 = [=] _sp_hybrid (const index_type& iii) mutable
        {
            auto i = iii;
            int iv = i.lb();
            i.lb() = the_lb_glob;
            
            ctrs::array<int, 3> offst(0,0,0);
            grid::node_idx_t inode = grid::cell_to_node(i, offst);
            std::size_t voffst = 3*(inode.i() + nni*inode.j() + nni*nnj*inode.k());
            const auto xx = g_img.get_coords(inode);
            for (int d = 0; d < 3; ++d) cptr[voffst+d] = xx[d];
        };
        dispatch::execute(range,  load0);
        dispatch::execute(range2, load1);
        
        char* data_base  = (char*)(&vdr[0]);
        char* coord_base = (char*)(&cdr[0]);
        
        std::size_t data_size  = data_raw.size()*sizeof(data_t);
        std::size_t coord_size = coord_raw.size()*sizeof(coor_t);
        
        //In characters
        const auto b64_size = [](std::size_t sz)
        {
            std::size_t base_size = 4*(sz/3);
            std::size_t extra     = sz % 3;
            if (extra > 0) base_size += 4;
            return base_size;
        };
        
        //In characters
        const auto b64_bulk_range = [](std::size_t sz)
        {
            std::size_t base_size = (sz/3);
            return base_size;
        };
        
        //In characters
        const auto b64_resid = [](std::size_t sz)
        {
            std::size_t base_size = sz%3;
            return base_size;
        };
        
        std::size_t header_size = b64_size(sizeof(unsigned int));
        std::size_t var_size_b64 = header_size + b64_size(data_size/num_vars);
        std::size_t cor_size_b64 = header_size + b64_size(coord_size);
        
        std::size_t var_bulk = b64_bulk_range(data_size/num_vars);
        std::size_t crd_bulk = b64_bulk_range(coord_size);
        
        std::size_t var_resid = b64_resid(data_size/num_vars);
        std::size_t crd_resid = b64_resid(coord_size);
        
        const auto bulk_range = dispatch::ranges::make_range(0UL, crd_bulk + num_vars*var_bulk);
        
        std::size_t var_bytes = ni*nj*nk*sizeof(data_t);
        
        auto target_img = utils::make_vec_image(target);
        auto load2 = [=] _sp_hybrid (const std::size_t& idx) mutable
        {
            constexpr char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
            constexpr char nil = '=';
            char* read_addr = coord_base + 3*idx;
            std::size_t o_idx = header_size + 4*idx;
            bool iscoord = true;
            std::size_t index = 0;
            std::size_t var = 0;
            if (idx >= crd_bulk)
            {
                std::size_t delta = idx - crd_bulk;
                var   = delta/var_bulk;
                index = delta % var_bulk;
                o_idx = cor_size_b64 + var*var_size_b64 + 4*index + header_size;
                read_addr =  data_base;
                read_addr += var*var_bytes;
                read_addr += 3*index;
                iscoord = false;
            }
            
            // if (!iscoord && ( >= data_size))
            // {
                
            // }
            
            char b0, b1, b2, c0, c1, c2, c3;
            b0 = read_addr[0];
            b1 = read_addr[1];
            b2 = read_addr[2];
            c0 = ((b0 & 0b11111100) >> 2);
            c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
            c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
            c3 = ((b2 & 0b00111111));
            
            target_img[o_idx + 0] = table[c0];
            target_img[o_idx + 1] = table[c1];
            target_img[o_idx + 2] = table[c2];
            target_img[o_idx + 3] = table[c3];
            
            std::size_t resid = iscoord?crd_resid:var_resid;
            if (resid > 0)
            {
                char* read_end        = coord_base + 3*crd_bulk;
                std::size_t write_idx = header_size + 4*crd_bulk;
                if (!iscoord)
                {
                    read_end = data_base + var*var_bytes + 3*var_bulk;
                    write_idx = cor_size_b64 + var*var_size_b64 + header_size + 4*var_bulk;
                }
                char last3[3]    = {0, 0, 0};
                char last4[4]    = {0, 0, 0, 0};
                char writebuf[4] = {nil, nil, nil, nil};
                int num2write[3] = {0, 2, 3};
                int numTerms[3] = {0, 2, 1};
                for (int i = 0; i < resid; i++)
                {
                    last3[i] = read_end[i];
                }
                b0 = last3[0];
                b1 = last3[1];
                b2 = last3[2];
                c0 = ((b0 & 0b11111100) >> 2);
                c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
                c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
                c3 = ((b2 & 0b00111111));
                last4[0] = c0;
                last4[1] = c1;
                last4[2] = c2;
                last4[3] = c3;
                for (int i = 0; i < num2write[resid]; ++i)
                {
                    writebuf[i] = table[last4[i]];
                }
                
                target_img[write_idx + 0] = writebuf[0];
                target_img[write_idx + 1] = writebuf[1];
                target_img[write_idx + 2] = writebuf[2];
                target_img[write_idx + 3] = writebuf[3];
            }
        };
        
        dispatch::execute(bulk_range, load2, arr.device());
        obuf.itransfer();
        
        std::stringstream ssv, ssc;
        unsigned int vbytes = data_raw.size()  * sizeof(data_t) / num_vars;
        unsigned int cbytes = coord_raw.size() * sizeof(coor_t);
        stream_base_64(ssv, &vbytes, 1);
        stream_base_64(ssc, &cbytes, 1);
        std::string sc = ssc.str();
        std::string sv = ssv.str();
        
        //Add the headers
        for (int i = 0; i < header_size; ++i)
        {
            obuf[i] = sc[i];
            for (std::size_t v = 0; v < num_vars; ++v)
            {
                std::size_t offst = v*var_size_b64 + cor_size_b64;
                obuf[offst + i] = sv[i];
            }
        }
        
        /*
        std::vector<float> allcrd = coord_raw.data(arr.device());
        std::stringstream ss;
        stream_base_64(ss, allcrd);
        std::ofstream f("crd" + std::to_string(the_lb));
        f << ss.str();
        */
        
        return std::make_pair(var_size_b64, cor_size_b64);
    }
}