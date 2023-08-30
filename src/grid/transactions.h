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
    };
    
    template <const std::size_t gdim>
    struct patch_fill_t
    {
        // for construction of this object, see get_transaction.h
        constexpr static int max_size = static_math::pow<2, gdim>::value;
        
        grid_rect_copy_t patches;
        ctrs::array<int, 3> i_skip;
        ctrs::array<ctrs::array<int, 3>, max_size> delta_i;
        int num_oslot;
        
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
            ctrs::array<element_type, max_size> oslot;
            constexpr std::size_t elem_size = sizeof(element_type);
            element_type elem;
            char* raw_e_addr = (char*)&oslot;
            std::size_t output = 0;
            grid::cell_idx_t idx;
            idx.lb() = patches.source.min(3);
            for (idx.k() = patches.source.min(2); idx.k() < patches.source.max(2); idx.k() += i_skip[2]){
            for (idx.j() = patches.source.min(1); idx.j() < patches.source.max(1); idx.j() += i_skip[1]){
            for (idx.i() = patches.source.min(0); idx.i() < patches.source.max(0); idx.i() += i_skip[0]){
                elem = f_val_t(0.0);
                algs::static_for<0, max_size>([&](const auto& ii)
                {
                    constexpr int nn = ii.value;
                    auto donor = idx;
                    donor.i() += delta_i[nn][0];
                    donor.j() += delta_i[nn][1];
                    donor.k() += delta_i[nn][2];
                    elem += array.get_elem(donor);
                });
                
                for (auto& k: oslot) k = coeff*elem;
                std::copy(raw_e_addr, raw_e_addr + num_oslot*elem_size, buf + output);
                output += num_oslot*elem_size;
            }}}
            return output;
        }

        template <typename arr_t>
        _finline_ std::size_t extract(arr_t& array, const char* buf) const
        {
            return patches.extract(array, buf);
        }
    };
}