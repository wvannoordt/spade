#pragma once

#include "core/bounding_box.h"
#include "core/static_math.h"

namespace spade::grid
{
    //Note: this will take some significant re-working for GPU implementation!

    //transaction pattern for a rectangular copy of a region of cells,
    //used for direct injection for cartesian configurations or
    //same-level transactions on AMR grids
    struct grid_rect_copy_t
    {
        bound_box_t<int, 4> source, dest;
        int rank_recv, rank_send;
        std::size_t send_volume() const { return source.volume(); }
        std::size_t recv_volume() const { return dest.volume(); }
        
        int receiver() const { return rank_recv; }
        int sender()   const { return rank_send; }
        
        template <typename arr_t>
        _finline_ std::size_t insert(const arr_t& array, char* buf) const
        {
            using element_type = typename arr_t::alias_type;
            constexpr std::size_t elem_size = sizeof(element_type);
            element_type elem;
            char* raw_e_addr = (char*)&elem;
            std::size_t output = 0;
            grid::cell_idx_t idx;
            algs::md_loop(idx, source, [&](const auto& icell)
            {
                elem = array.get_elem(idx);
                std::copy(raw_e_addr, raw_e_addr + elem_size, buf + output);
                output += elem_size;
            });
            return output;
        }
        
        template <typename arr_t>
        _finline_ std::size_t extract(arr_t& array, const char* buf) const
        {
            using element_type = typename arr_t::alias_type;
            element_type elem;
            char* raw_e_addr = (char*)&elem;
            std::size_t elem_size = sizeof(element_type);
            std::size_t output = 0;
            grid::cell_idx_t idx;
            algs::md_loop(idx, dest, [&](const auto& icell)
            {
                std::copy(buf + output, buf+output+elem_size, raw_e_addr);
                array.set_elem(idx, elem);
                output += elem_size;
            });
            return output;
        }
        
        bool reducible() const { return false; }
        auto reduce() const { return *this; }
    };
    
    template <const std::size_t gdim>
    struct patch_fill_t
    {
        // for construction of this object, see get_transaction.h
        constexpr static int max_size = static_math::pow<2, gdim>::value;
        
        grid_rect_copy_t patches;
        ctrs::array<int, 3> i_coeff, i_incr;
        
        //i_incr[d] = 1 if donor is finer in direction d
        //i_incr[d] = 0 otherwise
        
        //i_coeff[d] = 1 if donor is coarser in direction d
        //i_coeff[d] = 0 otherwise
        
        //NOTE: this type of exchange will "pack" the sending patch
        //into memory of the size expected to receive, hence the anomaly
        //in send_volume()
        std::size_t send_volume() const { return patches.recv_volume(); }
        std::size_t recv_volume() const { return patches.recv_volume(); }
        
        int receiver() const { return patches.receiver(); }
        int sender()   const { return patches.sender(); }
        
        template <typename arr_t>
        _finline_ std::size_t insert(const arr_t& array, char* buf) const
        {
            using element_type   = typename arr_t::alias_type;
            using f_val_t        = typename arr_t::fundamental_type;
            constexpr auto coeff = f_val_t(1.0)/max_size;
            const int lb = patches.source.min(3);
            constexpr std::size_t elem_size = sizeof(element_type);
            element_type elem = 0.0;
            char* rptr = (char*)(&elem);
            
            std::size_t output = 0;
            for (int dk = 0; dk < patches.dest.size(2); ++dk){
            for (int dj = 0; dj < patches.dest.size(1); ++dj){
            for (int di = 0; di < patches.dest.size(0); ++di){
                //each cell in the output range
                elem = 0.0;
                // grid::cell_idx_t recvr;
                // recvr.lb() = patches.dest.min(3);
                // recvr.i() = patches.dest.min(0) + di;
                // recvr.j() = patches.dest.min(1) + dj;
                // recvr.k() = patches.dest.min(2) + dk;
                
                
                // using t_t = spade::grid::cell_idx_t;
                // bool debug = (recvr == t_t(-1, 0, 0, 40));
                // if (debug) print("DEBUG");
                algs::static_for<0,max_size>([&](const auto iii)
                {
                    constexpr int iv = iii.value;
                    constexpr int d0 = (iv >> 0)&1;
                    constexpr int d1 = (iv >> 1)&1;
                    constexpr int d2 = (iv >> 2)&1;
                    
                    int i0 = di<<(i_coeff[0]+1);
                    int j0 = dj<<(i_coeff[1]+1);
                    int k0 = dk<<(i_coeff[2]+1);
                    i0 = i0 >> 1;
                    j0 = j0 >> 1;
                    k0 = k0 >> 1;
                    
                    
                    const int donor_di = i0 + d0*i_incr[0];
                    const int donor_dj = j0 + d1*i_incr[1];
                    const int donor_dk = k0 + d2*i_incr[2];
                    
                    grid::cell_idx_t idx;
                    idx.lb() = lb;
                    idx.i()  = patches.source.min(0) + donor_di;
                    idx.j()  = patches.source.min(1) + donor_dj;
                    idx.k()  = patches.source.min(2) + donor_dk;
                    
                    // if (debug) print(idx);
                    
                    elem += array.get_elem(idx);
                });
                // if (debug) print("i_incr",  i_incr);
                // if (debug) print("i_coeff", i_coeff);
                // if (debug) print("DONE");
                
                
                elem *= coeff;
                std::copy(rptr, rptr+elem_size, buf + output);
                output += elem_size;
            }}}
            return output;
        }

        template <typename arr_t>
        _finline_ std::size_t extract(arr_t& array, const char* buf) const
        {
            std::size_t output = patches.extract(array, buf);
            return output;
        }
        
        bool reducible() const
        {
            bool output = true;
            for (const auto p: i_coeff)
            {
                output = output && (p == 0);
            }
            for (const auto p: i_incr)
            {
                output = output && (p == 0);
            } 
            return output;
        }
        auto reduce() const { return patches; }
        
    };
}