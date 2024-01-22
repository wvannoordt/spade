#pragma once

#include "grid/grid_index_types.h"

#include "omni/info.h"
#include "omni/stencil.h"
#include "omni/info_union.h"
#include "omni/stencil_union.h"

namespace spade::omni
{
    // Wraps a stencil-data instance with a type that 
    // effectively re-numbers the elements for the sake of 
    // nesting.
    // Note that input_t could be an alias itself!
    template <typename interpret_t, typename input_t>
    struct stencil_alias_t
    {
        using stencil_type        = interpret_t; // used for any further (deeper) aliases
        using native_stencil_type = typename input_t::stencil_type;
        const input_t& native;
        _sp_hybrid stencil_alias_t(const input_t& native_in) : native{native_in}{}

        _sp_hybrid auto& root() {return native.root();}
        _sp_hybrid const auto& root() const {return native.root();}

        //const qualified
        template <const grid::array_centering ctr, udci::integral_t ii>
        _sp_hybrid decltype(auto) seek_element(const udci::idx_const_t<ii>& idx) const
        {
            // strategy: generate as much heat as possible during compilation so the
            // data center doesn't get too cold
            using offset_type     = offset_at<interpret_t, ctr, ii>;
            constexpr int trans_i = index_of<native_stencil_type, offset_type>;
            return native.template seek_element<ctr>(udci::idx_const_t<trans_i>());
        }
        
        //not const qualified
        template <const grid::array_centering ctr, udci::integral_t ii>
        _sp_hybrid decltype(auto) seek_element(const udci::idx_const_t<ii>& idx)
        {
            using offset_type     = offset_at<interpret_t, ctr, ii>;
            constexpr int trans_i = index_of<native_stencil_type, offset_type>;
            return native.template seek_element<ctr>(udci::idx_const_t<trans_i>());
        }
        
        //const qualified
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_cell())
        _sp_hybrid constexpr decltype(auto) cell(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::cell_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_face())
        _sp_hybrid constexpr decltype(auto) face(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::face_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_node())
        _sp_hybrid constexpr decltype(auto) node(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::node_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_edge())
        _sp_hybrid constexpr decltype(auto) edge(const udci::idx_const_t<ii>& idx) const
        {
            return seek_element<grid::edge_centered>(idx);
        }
        
        //not const qualified
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_cell())
        _sp_hybrid constexpr decltype(auto) cell(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::cell_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_face())
        _sp_hybrid constexpr decltype(auto) face(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::face_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_node())
        _sp_hybrid constexpr decltype(auto) node(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::node_centered>(idx);
        }
        
        template <udci::integral_t ii>
        requires (ii < interpret_t::num_edge())
        _sp_hybrid constexpr decltype(auto) edge(const udci::idx_const_t<ii>& idx)
        {
            return seek_element<grid::edge_centered>(idx);
        }

    };

    template <typename interpret_t, typename input_t>
    _sp_hybrid auto interpret_stencil(const input_t& input)
    {
        return stencil_alias_t<interpret_t, input_t>(input);
    }

    template <typename interpret_t, typename input_t>
    requires (std::same_as<typename input_t::stencil_type, interpret_t>)
    _sp_hybrid auto interpret_stencil(const input_t& input)
    {
        return input;
    }
}