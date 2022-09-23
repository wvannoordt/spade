#pragma once

#include <tuple>
#include <type_traits>

#include "core/grid.h"
#include "core/coord_system.h"

namespace spade::fetch
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
    
    template <typename data_t> struct face_state
    {
        data_t data;
    };
    
    template <typename data_t> struct face_state_grad
    {
        data_t data;
    };
    
    template <typename data_t> struct face_normal
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
        const static std::size_t num_params = sizeof...(infos_t);
        std::tuple<infos_t...> elements;
        template <const std::size_t idx> const auto& item(void) const {return std::get<idx>(elements).data;}
        template <const std::size_t idx> auto& item(void) {return std::get<idx>(elements).data;}
    };
    
    template <typename cell_info_t> struct left_right
    {
        cell_info_t left, right;
    };
    
    template <const std::size_t stencil_size, typename cell_info_t> struct flux_line
    {
        ctrs::array<cell_info_t, stencil_size> stencil;
    };
    
    template <const std::size_t stencil_size, typename cell_info_t>
    static std::ostream & operator<<(std::ostream & os, const flux_line<stencil_size, cell_info_t>& ftch)
    {
        os << "Line-stencil, size " << stencil_size <<".\nData:\n";
        for (auto i: range(stencil_size))
        {
            os << ftch[i] << "\n";
        }
        return os;
    }
    
    template <typename cell_stencil_t, typename face_info_t> struct face_fetch_t
    {
        cell_stencil_t cell_data;
        face_info_t    face_data;
    };
    
    template <typename cell_stencil_t, typename face_info_t>
    static std::ostream & operator<<(std::ostream & os, const face_fetch_t<cell_stencil_t, face_info_t>& ftch)
    {
        os << "Face fetch:\n >> Cell data:\n" << ftch.cell_data << "\n >> Face data:\n" << ftch.face_data;
        return os;
    }
    
    namespace detail
    {
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_single_cell_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const grid::cell_idx_t& icell,
            const int& idir,
            cell_state<output_t>& output)
        {
            for (auto n: range(0,output.data.size())) output.data[n] = prims(n, icell[0], icell[1], icell[2], icell[3]);
        }
        
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_single_cell_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const grid::cell_idx_t& icell,
            const int& idir,
            cell_normal<output_t>& output)
        {
            const ctrs::array<typename grid_t::dtype,3> xyz_c = ar_grid.get_comp_coords(icell);
            output.data = coords::calc_normal_vector(ar_grid.coord_sys(), xyz_c, icell, idir);
        }
        
        template <typename cell_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_single_cell_data(
            const grid_t& ar_grid,
            const array_t& prims,
            const grid::cell_idx_t& icell,
            const int& idir,
            cell_info_t& info)
        {
            auto& p = info.elements;
            static_for<0,cell_info_t::num_params>([&](auto i) -> void
            {
                get_single_cell_info_value(ar_grid, prims, icell, idir, std::get<i.value>(info.elements));
            });
        }
        
        template <typename cell_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_cell_stencil_data(const grid_t& ar_grid, const array_t& prims, const grid::face_idx_t& iface, left_right<cell_info_t>& lr)
        {
            const int idir = grid::get_face_dir(iface);
            grid::cell_idx_t il = grid::face_to_cell(iface, 0);
            grid::cell_idx_t ir = grid::face_to_cell(iface, 1);
            get_single_cell_data(ar_grid, prims, il, idir, lr.left);
            get_single_cell_data(ar_grid, prims, ir, idir, lr.right);
        }
        
        template <const std::size_t flux_size, typename cell_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_cell_stencil_data(const grid_t& ar_grid, const array_t& prims, const grid::face_idx_t& iface, flux_line<flux_size, cell_info_t>& sten)
        {
            const int idir = grid::get_face_dir(iface);
            static_for<0,flux_size>([&](auto i) -> void
            {
                grid::cell_idx_t icell = grid::face_to_cell(iface, 0);
                icell[idir] += i.value - (int)(flux_size/2) + 1;
                get_single_cell_data(ar_grid, prims, icell, idir, sten.stencil[i.value]);
            });
        }
        
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_face_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const grid::face_idx_t& iface,
            face_normal<output_t>& output)
        {
            const ctrs::array<typename grid_t::dtype,3> xyz_c = ar_grid.get_comp_coords(iface);
            output.data = coords::calc_normal_vector(ar_grid.coord_sys(), xyz_c, iface, grid::get_face_dir(iface));
        }
        
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_face_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const grid::face_idx_t& iface,
            face_state<output_t>& output)
        {
            const int idir = grid::get_face_dir(iface);
            const grid::cell_idx_t il = grid::face_to_cell(iface, 0);
            const grid::cell_idx_t ir = grid::face_to_cell(iface, 1);
            cell_state<output_t> ql, qr;
            get_single_cell_info_value(ar_grid, prims, il, idir, ql);
            get_single_cell_info_value(ar_grid, prims, ir, idir, qr);
            for (int n = 0; n < output.data.size(); ++n) output.data[n] = 0.5*ql.data[n] + 0.5*qr.data[n];
        }
        
        template <grid::multiblock_grid grid_t, ctrs::basic_array vec_t, typename idx_t>
        requires std::is_same<typename grid_t::coord_sys_type, coords::identity<typename grid_t::coord_type>>::value
        void transform_gradient(const grid_t& ar_grid, const idx_t& idx, ctrs::array<vec_t, 3>& data)
        {
            // do nothing (identity coordinates)
        }
        
        template <grid::multiblock_grid grid_t, ctrs::basic_array vec_t, typename idx_t>
        requires coords::diagonal_coordinate_system<typename grid_t::coord_sys_type>
        void transform_gradient(const grid_t& ar_grid, const idx_t& idx, ctrs::array<vec_t, 3>& data)
        {
            const auto& x   = ar_grid.get_comp_coords(idx);
            const auto& sys = ar_grid.coord_sys();
            const auto xi_x = 1.0/sys.xcoord.coord_deriv(x[0]);
            const auto xi_y = 1.0/sys.ycoord.coord_deriv(x[1]);
            const auto xi_z = 1.0/sys.zcoord.coord_deriv(x[2]);
            data[0] *= xi_x;
            data[1] *= xi_y;
            data[2] *= xi_z;
        }
        
        template <grid::multiblock_grid grid_t, ctrs::basic_array vec_t, typename idx_t>
        requires coords::dense_coordinate_system<typename grid_t::coord_sys_type>
        void transform_gradient(const grid_t& ar_grid, const idx_t& idx, ctrs::array<vec_t, 3> data)
        {
            const auto& x   = ar_grid.get_comp_coords(idx);
            const auto& sys = ar_grid.coord_sys();
            //invert 3x3 matrix
            auto deriv = sys.coord_deriv(x);
            linear_algebra::brute_solve_inplace(deriv, data);
            print("This isn't implemented yet", __FILE__, __LINE__);
            abort();
        }
        
        template <typename output_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_face_info_value(
            const grid_t& ar_grid,
            const array_t& prims,
            const grid::face_idx_t& iface,
            face_state_grad<output_t>& output)
        {
            const int idir0 = grid::get_face_dir(iface);
            const ctrs::array<int,3> idir(idir0, (idir0+1)%ar_grid.dim(), (idir0+2)%ar_grid.dim());
            const ctrs::array<typename grid_t::coord_type, 3> invdx
            (
                1.0/ar_grid.get_dx(0),
                1.0/ar_grid.get_dx(1),
                1.0/ar_grid.get_dx(2)
            );
            grid::cell_idx_t ic = grid::face_to_cell(iface, 0);
            
            output.data = 0.0;
            cell_state<typename output_t::value_type> q;
            auto apply_coeff_at = [&](const int& iset, const typename array_t::value_type& coeff, const grid::cell_idx_t& idx)
            {
                get_single_cell_info_value(ar_grid, prims, idx, idir[0], q);
                for (std::size_t i = 0; i < output.data[iset].size(); ++i) output.data[iset][i] += coeff*q.data[i]*invdx[iset];
            };
            
            apply_coeff_at(idir[0],  -1.0, ic);
            for (int ii = 1; ii < grid_t::dim(); ++ii)
            {
                ic[idir[ii]] += 1;
                apply_coeff_at(idir[ii],  0.25, ic);
                ic[idir[ii]] -= 2;
                apply_coeff_at(idir[ii], -0.25, ic);
                ic[idir[ii]] += 1;
            }
            ic[idir[0]] += 1;
            apply_coeff_at(idir[0],   1.0, ic);
            for (int ii = 1; ii < grid_t::dim(); ++ii)
            {
                ic[idir[ii]] += 1;
                apply_coeff_at(idir[ii],  0.25, ic);
                ic[idir[ii]] -= 2;
                apply_coeff_at(idir[ii], -0.25, ic);
                ic[idir[ii]] += 1;
            }
            
            const auto& x_face = ar_grid.get_comp_coords(iface);
            transform_gradient(ar_grid, iface, output.data);
        }
        
        
        template <typename face_info_t, grid::multiblock_grid grid_t, grid::multiblock_array array_t>
        void get_face_data(const grid_t& grid, const array_t& prims, const grid::face_idx_t& iface, face_info_t& face_data)
        {
            auto& p = face_data.elements;
            static_for<0,face_info_t::num_params>([&](auto i) -> void
            {
                get_face_info_value(grid, prims, iface, std::get<i.value>(p));
            });
        }
    }
    
    template <grid::multiblock_grid grid_t, grid::multiblock_array array_t, typename flux_in_t>
    void get_flux_data(const grid_t& grid, const array_t& prims, const grid::face_idx_t& iface, flux_in_t& flux_input)
    {
        detail::get_cell_stencil_data(grid, prims, iface, flux_input.cell_data);
        detail::get_face_data(grid, prims, iface, flux_input.face_data);
    }
}