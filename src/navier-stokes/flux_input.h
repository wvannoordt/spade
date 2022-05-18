#pragma once

#include <tuple>

namespace cvdf::flux_input
{
    template <typename data_t> struct cell_state
    {
        data_t data;
    };
    
    template <typename data_t> struct cell_normal
    {
        data_t data;
    };
    
    template <typename data_t> struct cell_index
    {
        data_t data;
    };
    
    template <typename data_t> struct cell_coord
    {
        data_t data;
    };
    
    template <typename... infos_t> struct cell_info
    {
        const static std::size_t num_params = sizeof...(infos_t);
        std::tuple<infos_t...> elements;
        template <const std::size_t idx> const auto& item(void) const {return std::get<idx>(elements).data;}
        template <const std::size_t idx> auto& item(void) {return std::get<idx>(elements).data;}
    };
    
    template <typename... infos_t> struct face_info
    {
        std::tuple<infos_t...> data;
    };
    
    template <typename cell_info_t> struct left_right
    {
        cell_info_t left, right;
    };
    
    template <const std::size_t stencil_size, typename cell_info_t> struct flux_line
    {
        ctrs::array<cell_info_t, stencil_size> stencil;
    };
    
    template <typename cell_stencil_t, typename face_info_t> struct flux_input_t
    {
        cell_stencil_t cell_data;
        face_info_t    face_data;
    };
    
    namespace detail
    {
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_single_cell_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const ctrs::array<grid::cell_t<int>, 4>& icell,
            const int& idir,
            cell_state<output_t>& output)
        {
            for (auto n: range(0,output.data.size())) output.data[n[0]] = prims(n[0], icell[0], icell[1], icell[2], icell[3]);
        }
        
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_single_cell_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const ctrs::array<grid::cell_t<int>, 4>& icell,
            const int& idir,
            cell_normal<output_t>& output)
        {
            const ctrs::array<typename grid_t::dtype,3> xyz_c = ar_grid.get_comp_coords(icell[0], icell[1], icell[2], icell[3]);
            output.data = coords::calc_normal_vector(ar_grid.coord_sys(), xyz_c, icell, idir);
        }
        
        template <typename cell_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_single_cell_data(
            const grid_t& ar_grid,
            const array_t& prims,
            const ctrs::array<grid::cell_t<int>, 4>& icell,
            const int& idir,
            cell_info_t& info)
        {
            auto& p = info.elements;
            // static_for<0,std::tuple_size<decltype(info.elements)>([&](auto i) -> void
            static_for<0,cell_info_t::num_params>([&](auto i) -> void
            {
                get_single_cell_info_value(ar_grid, prims, icell, idir, std::get<i.value>(info.elements));
            });
        }
        
        template <typename cell_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_cell_stencil_data(const grid_t& ar_grid, const array_t& prims, const ctrs::array<grid::face_t<int>, 5>& iface, left_right<cell_info_t>& lr)
        {
            const int idir = iface[0];
            ctrs::array<grid::cell_t<int>, 4> il((int)iface[1], (int)iface[2], (int)iface[3], (int)iface[4]);
            ctrs::array<grid::cell_t<int>, 4> ir((int)iface[1], (int)iface[2], (int)iface[3], (int)iface[4]);
            ir[idir] += 1;
            get_single_cell_data(ar_grid, prims, il, idir, lr.left);
            get_single_cell_data(ar_grid, prims, ir, idir, lr.right);
        }
        
        template <const std::size_t flux_size, typename cell_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_cell_stencil_data(const grid_t& ar_grid, const array_t& prims, const ctrs::array<grid::face_t<int>, 5>& iface, flux_line<flux_size, cell_info_t>& sten)
        {
            const int idir = iface[0];
            static_for<0,flux_size>([&](auto i) -> void
            {
                ctrs::array<grid::cell_t<int>, 4> icell((int)iface[1], (int)iface[2], (int)iface[3], (int)iface[4]);
                icell[idir] += i.value - (int)(flux_size/2) + 1;
                get_single_cell_data(ar_grid, prims, icell, idir, sten.stencil[i.value]);
            });
        }
        
        template <typename face_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_face_stencil_data(const grid_t& grid, const array_t& prims, const ctrs::array<grid::face_t<int>, 5>& iface, face_info_t& face_data)
        {
            //TODO [viscous implementation]
        }
    }
    
    template <grid::multiblock_grid grid_t, grid::multiblock_array array_t, typename flux_in_t>
    void get_flux_data(const grid_t& grid, const array_t& prims, const ctrs::array<grid::face_t<int>, 5>& iface, flux_in_t& flux_input)
    {
        detail::get_cell_stencil_data(grid, prims, iface, flux_input.cell_data);
        detail::get_face_stencil_data(grid, prims, iface, flux_input.face_data);
    }
}