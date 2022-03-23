#pragma once

#include <iostream>
#include <fstream>
#include <filesystem>

#include "grid.h"
#include "print.h"
#include "utils.h"

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
            out_str << ((cvdf_dim==3)*(1)+obj.get_num_blocks(2)*obj.get_num_cells(2)) << "\n";
            auto bnd = obj.get_bounds();
            out_str << "ORIGIN "  << bnd.min(0)    << " " << bnd.min(1)    << " " << bnd.min(2)    << "\n";
            out_str << "SPACING " << obj.get_dx(0) << " " << obj.get_dx(1) << " " << (cvdf_dim==3)*obj.get_dx(2) << "\n";
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
            std::size_t num_node_centered_arrays = utils::reduce_over_params(test_node_center, arrays...);
            
            auto test_cell_center = [](auto i) -> std::size_t {return (i.centering_type()==grid::cell_centered)?1:0;};
            std::size_t num_cell_centered_arrays = utils::reduce_over_params(test_cell_center, arrays...);
            
            int ct = 0;
            auto current_centering = grid::node_centered;
            auto output_array_data = [&](auto arr) -> void
            {
                if (arr.centering_type()==current_centering)
                {
                    for (auto j: range(0,arr.get_major_dims().total_size()))
                    {
                        for (auto i: range(0,arr.get_minor_dims().total_size()))
                        {
                            std::string name = "var" + zfill(ct++, 5);
                            out_str << "SCALARS " << name << " double\nLOOKUP_TABLE default\n";
                            auto grid_range = obj.get_range(arr.centering_type(), grid::exclude_exchanges);
                            auto r1 = grid_range.subrange(0)*range(0, obj.get_num_blocks(0));
                            auto r2 = grid_range.subrange(1)*range(0, obj.get_num_blocks(1));
                            auto r3 = grid_range.subrange(2)*range(0, obj.get_num_blocks(2));
                            auto total_range = r1*r2*r3;
                            
                            for (auto midx: total_range)
                            {
                                std::size_t lb = obj.collapse_block_num(midx[1], midx[3], midx[5]);
                                out_str << arr.unwrap_idx(i[0], midx[0], midx[2], midx[4], lb, j[0]) << "\n";
                            }
                        }
                    }
                }
            };
            
            if (num_node_centered_arrays>0) out_str << "POINT_DATA " << obj.get_range(grid::node_centered, grid::exclude_exchanges).size() << "\n";
            utils::foreach_param(output_array_data, arrays...);
            current_centering = grid::cell_centered;
            if (num_cell_centered_arrays>0) out_str << "CELL_DATA "  << obj.get_range(grid::cell_centered, grid::exclude_exchanges).size() << "\n";
            utils::foreach_param(output_array_data, arrays...);
        }
        
        std::string ntab(const std::size_t& n, const std::size_t& tab_size = 4) { return std::string(tab_size*n, ' '); }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t>
        void output_parallel_header_file(
            output_stream_t& out_str,
            const std::string& block_proto_string,
            const std::size_t& num_zeros,
            const grid_output_t& obj)
        {
            out_str << "<?xml version=\"1.0\"?>\n";
            out_str << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n";
            out_str << ntab(1) << "<vtkMultiBlockDataSet>\n";
            out_str << ntab(2) << "<Block index =\"0\">\n";
            // <DataSet index="100" file="flow/lev000/flow_bk0000100.vtr"/>
            for (auto i: range(0, obj.get_num_global_blocks()))
            {
                out_str << ntab(3) << utils::strformat("<DataSet index=\"{}\" file=\"{}\"/>\n", i[0], utils::strformat(block_proto_string, utils::zfill(i[0],num_zeros)));
            }
            out_str << ntab(2) << "</Block>\n";
            out_str << ntab(1) << "</vtkMultiBlockDataSet>\n";
            out_str << "</VTKFile>\n";
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, typename... arrays_t>
        void output_parralel_block_file(output_stream_t& out_str, const grid_output_t& obj, arrays_t... arrays)
        {
            out_str << "hello\n";
        }
        
        template <grid::multiblock_grid grid_output_t, typename... arrays_t>
        std::string output_grid_parallel(const std::string& out_dir, const std::string& out_name_no_extension, const grid_output_t& obj, arrays_t... arrays)
        {
            std::filesystem::path out_path(out_dir);
            if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
            std::filesystem::path header_filename(out_dir);
            header_filename /= (out_name_no_extension + ".vtm");
            std::size_t numzeros = 7;
            
            std::filesystem::path block_pth(out_dir);
            block_pth /= ("data_" + out_name_no_extension);
            if (!std::filesystem::is_directory(block_pth)) std::filesystem::create_directory(block_pth);
            block_pth /= "block{}.vtr";
            std::string block_template = block_pth;
            
            if (obj.group().isroot())
            {
                std::ofstream header_file_strm(header_filename);
                
                output_parallel_header_file(header_file_strm, block_template, numzeros, obj);
            }
            for (auto i: range(0, obj.get_num_local_blocks()))
            {
                std::string block_file_name = utils::strformat(block_template, zfill(obj.get_partition().get_global_block(i[0]), numzeros));
                std::ofstream block_file_strm(block_file_name);
                output_parralel_block_file(block_file_strm, obj, arrays...);
            }
            return header_filename;
        }
    }
    
    template <grid::multiblock_grid grid_output_t, typename... arrays_t>
    std::string output_vtk(const std::string& out_dir, const std::string& out_name_no_extension, const grid_output_t& obj, arrays_t... arrays)
    {
        if (obj.group().size()>1)
        {
            return detail::output_grid_parallel(out_dir, out_name_no_extension, obj, arrays...);
        }
        else
        {
            std::filesystem::path out_path(out_dir);
            if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
            out_path /= (out_name_no_extension + ".vtk");
            std::string total_filename = out_path;
            std::ofstream out_str(total_filename);
            detail::output_grid_serial(out_str, obj, arrays...);
            return total_filename;
        }
    }
}