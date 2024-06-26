#pragma once

#include "core/bounding_box.h"
#include "core/static_math.h"

namespace spade::grid
{
    
    enum transaction_tag_t
    {
        null_transaction       = 0,
        zface_transaction      = 1,
        yface_transaction      = 2,
        xface_transaction      = 3,
        zhalf_face_transaction = 4,
        yhalf_face_transaction = 5,
        xhalf_face_transaction = 6,
        zqutr_face_transaction = 7,
        yqutr_face_transaction = 8,
        xqutr_face_transaction = 9,
        zedge_transaction      = 10,
        yedge_transaction      = 11,
        xedge_transaction      = 12,
        zhalf_edge_transaction = 13,
        yhalf_edge_transaction = 14,
        xhalf_edge_transaction = 15,
        crnr_transaction       = 16
    };
    
    //To-do
    // static_assert(
    //     utils::is_decreasing<
    //         crnr_transaction,
    //         half_edge_transaction,
    //         edge_transaction,
    //         qutr_face_transaction,
    //         half_face_transaction,
    //         face_transaction,
    //         null_transaction
    //     >, "Invalid ordering for transaction varieties!");
    
    inline std::string trs2str(const transaction_tag_t& tg)
    {
        switch (tg)
        {
            case null_transaction:       return "null_transaction";
            case zface_transaction:      return "zface_transaction";
            case yface_transaction:      return "yface_transaction";
            case xface_transaction:      return "xface_transaction";
            case zhalf_face_transaction: return "zhalf_face_transaction";
            case yhalf_face_transaction: return "yhalf_face_transaction";
            case xhalf_face_transaction: return "xhalf_face_transaction";
            case zqutr_face_transaction: return "zqutr_face_transaction";
            case yqutr_face_transaction: return "yqutr_face_transaction";
            case xqutr_face_transaction: return "xqutr_face_transaction";
            case zedge_transaction:      return "zedge_transaction";
            case yedge_transaction:      return "yedge_transaction";
            case xedge_transaction:      return "xedge_transaction";
            case zhalf_edge_transaction: return "zhalf_edge_transaction";
            case yhalf_edge_transaction: return "yhalf_edge_transaction";
            case xhalf_edge_transaction: return "xhalf_edge_transaction";
            case crnr_transaction:       return "crnr_transaction";
        }
        return "trs2str ERR";
    }
    
    inline std::ostream& operator << (std::ostream& os, const transaction_tag_t& tg)
    {
        os << trs2str(tg);
        return os;
    }

    //transaction pattern for a rectangular copy of a region of cells,
    //used for direct injection for cartesian configurations or
    //same-level transactions on AMR grids
    struct grid_rect_copy_t
    {
        transaction_tag_t tag = null_transaction;
        bound_box_t<int, 4> source, dest;
        int rank_recv, rank_send;
        
        std::size_t glob_source_blk, glob_dest_blk;
        
        std::size_t glob_source_block() const { return glob_source_blk; }
        std::size_t glob_dest_block()   const { return glob_dest_blk; }
        
        std::size_t send_volume() const { return source.volume(); }
        std::size_t recv_volume() const { return dest.volume(); }
        
        int receiver() const { return rank_recv; }
        int sender()   const { return rank_send; }
        
        transaction_tag_t get_tag() const { return tag; }
        
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
        transaction_tag_t tag = null_transaction;
        
        // for construction of this object, see get_transaction.h
        constexpr static int max_size = static_math::pow<2, gdim>::value;
        
        grid_rect_copy_t patches;
        std::size_t glob_source_block() const { return patches.glob_source_blk; }
        std::size_t glob_dest_block()   const { return patches.glob_dest_blk; }
        
        ctrs::array<int, 3> i_coeff, i_incr;
        
        //NOTE: this type of exchange will "pack" the sending patch
        //into memory of the size expected to receive, hence the anomaly
        //in send_volume()
        std::size_t send_volume() const { return patches.recv_volume(); }
        std::size_t recv_volume() const { return patches.recv_volume(); }
        
        int receiver() const { return patches.receiver(); }
        int sender()   const { return patches.sender(); }
        
        transaction_tag_t get_tag() const { return tag; }
        
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
                    
                    elem += array.get_elem(idx);
                });
                
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