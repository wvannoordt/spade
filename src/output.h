#pragma once

#include <iostream>
#include <fstream>

#include "grid.h"
#include "print.h"
#include "algs.h"

namespace cvdf::output
{
    namespace detail
    {
        template <class output_stream_t> void output_serial_header(output_stream_t& out_str)
        {
            out_str << "# vtk DataFile Version 2.0\ncvdf output\nASCII\n";
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t>
        void output_mesh_data(output_stream_t& out_str, const grid_output_t& obj, const coords::identity<typename grid_output_t::dtype>& coord_sys)
        {
            out_str << "DATASET STRUCTURED_POINTS\nDIMENSIONS ";
            out_str << (1+obj.get_num_blocks(0)*obj.get_num_cells(0)) <<  " ";
            out_str << (1+obj.get_num_blocks(1)*obj.get_num_cells(1)) <<  " ";
            out_str << (1+obj.get_num_blocks(2)*obj.get_num_cells(2)) << "\n";
            auto bnd = obj.get_bounds();
            out_str << "ORIGIN "  << bnd.min(0)    << " " << bnd.min(1)    << " " << bnd.min(2)    << "\n";
            out_str << "SPACING " << obj.get_dx(0) << " " << obj.get_dx(1) << " " << obj.get_dx(2) << "\n";
            
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, coords::diagonal_coordinate_system diag_coord_sys_t>
        void output_mesh_data(output_stream_t& out_str, const grid_output_t& obj, const diag_coord_sys_t& coord_sys)
        {
            print("NOT IMPLEMENTED!", __FILE__, __LINE__);
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, typename... arrays_t>
        void output_grid_serial(output_stream_t& out_str, const grid_output_t& obj, arrays_t... arrays)
        {
            
            output_serial_header(out_str);
            output_mesh_data(out_str, obj, obj.coord_sys());
            
            auto test_node_center = [](auto i) -> std::size_t {return (i.centering_type()==grid::node_centered)?1:0;};
            std::size_t num_node_centered_arrays = algs::reduce_over_params(test_node_center, arrays...);
            
            auto test_cell_center = [](auto i) -> std::size_t {return (i.centering_type()==grid::cell_centered)?1:0;};
            std::size_t num_cell_centered_arrays = algs::reduce_over_params(test_cell_center, arrays...);
            
            int ct = 0;
            auto current_centering = grid::node_centered;
            auto output_array_data = [&](auto arr) -> void
            {
                if (arr.centering==current_centering)
                {
                    for (auto i: range(0,5))
                    {
                        std::string name = "var" + zfill(ct++, 5);
                        out_str << "SCALARS " << name << " double\nLOOKUP_TABLE default\n";
                        for (auto j: obj.get_range(arr.centering, grid::exclude_exchanges))
                        {
                            out_str << arr(i[0], j[0], j[1], j[2], j[3]) << "\n";
                        }
                    }
                }
            };
            
            if (num_node_centered_arrays>0) out_str << "POINT_DATA " << obj.get_range(grid::node_centered, grid::exclude_exchanges).size() << "\n";
            algs::foreach_param(output_array_data, arrays...);
            current_centering = grid::cell_centered;
            if (num_cell_centered_arrays>0) out_str << "CELL_DATA "  << obj.get_range(grid::cell_centered, grid::exclude_exchanges).size() << "\n";
            algs::foreach_param(output_array_data, arrays...);
        }
    }
    
    template <class output_stream_t, grid::multiblock_grid grid_output_t, typename... arrays_t>
    void output_grid(output_stream_t& out_str, const grid_output_t& obj, arrays_t... arrays)
    {
        if (obj.group().size()>1)
        {
            print("NOT IMPLEMENTED!", __FILE__, __LINE__);
        }
        else
        {
            detail::output_grid_serial(out_str, obj, arrays...);
        }
    }
}