#pragma once

#include <iostream>
#include <fstream>

#include "core/config.h"
#include "core/static_for.h"
#include "core/grid.h"
#include "core/print.h"
#include "core/utils.h"
#include "core/base_64.h"
#include "core/block_config.h"

namespace spade::io
{
    namespace detail
    {
        template <typename T> concept has_single_name    = requires(T t) {T::name();};
        template <typename T> concept has_multiple_names = requires(T t, const int& i) {T::name(i);};
        
        template <typename T> 
        requires (!has_single_name<T> && !has_multiple_names<T>)
        static std::string get_name(const int& i) {return "data";}
        
        template <typename T> 
        requires has_single_name<T>
        static std::string get_name(const int& i) {return T::name();}
        
        template <typename T> 
        requires has_multiple_names<T>
        static std::string get_name(const int& i) {return T::name(i);}
        
        
        template <typename arr_t> std::string get_array_name(const arr_t& arr, const int& i_minor, const int& i_major)
        {
            typedef typename arr_t::alias_type alias_t;
            bool has_provided_name = has_single_name<alias_t> || has_multiple_names<alias_t>;
            bool has_minor_suffix = !has_provided_name && arr.get_minor_dims().total_size()>1;
            bool has_major_suffix = arr.get_major_dims().total_size()>1;
            std::string base_name = "data";
            std::string minor_suffix = "";
            std::string major_suffix = "";
            if (has_minor_suffix) minor_suffix = utils::zfill(i_minor, 4);
            if (has_major_suffix) major_suffix = utils::zfill(i_major, 4);
            std::string join = "";
            if (has_provided_name) base_name = get_name<typename arr_t::alias_type>(i_minor);
            if (has_minor_suffix && has_major_suffix) join = "_";
            return base_name + minor_suffix + join + major_suffix;
        }
        
        template <class output_stream_t> static void output_serial_header(output_stream_t& out_str)
        {
            out_str << "# vtk DataFile Version 2.0\nspade output\nASCII\n";
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t>
        static void output_mesh_data(output_stream_t& out_str, const grid_output_t& obj, const coords::identity<typename grid_output_t::dtype>& coord_sys)
        {
            out_str << "DATASET STRUCTURED_POINTS\nDIMENSIONS ";
            out_str << (1+obj.get_num_blocks(0)*obj.get_num_cells(0)) <<  " ";
            out_str << (1+obj.get_num_blocks(1)*obj.get_num_cells(1)) <<  " ";
            out_str << ((obj.dim()==3)*(1)+obj.get_num_blocks(2)*obj.get_num_cells(2)) << "\n";
            auto bnd = obj.get_bounds();
            out_str << "ORIGIN "  << bnd.min(0)    << " " << bnd.min(1)    << " " << bnd.min(2)    << "\n";
            out_str << "SPACING " << obj.get_dx(0) << " " << obj.get_dx(1) << " " << (obj.dim()==3)*obj.get_dx(2) << "\n";
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, coords::diagonal_coordinate_system diag_coord_sys_t>
        static void output_mesh_data(output_stream_t& out_str, const grid_output_t& obj, const diag_coord_sys_t& coord_sys)
        {
            out_str << "DATASET RECTILINEAR_GRID\nDIMENSIONS ";
            int nx = (1+obj.get_num_blocks(0)*obj.get_num_cells(0));
            int ny = (1+obj.get_num_blocks(1)*obj.get_num_cells(1));
            int nz = ((obj.dim()==3)*(1)+obj.get_num_blocks(2)*obj.get_num_cells(2));
            out_str << nx <<  " " << ny <<  " " << nz << "\n";
            auto bnd = obj.get_bounds();
            out_str << "X_COORDINATES " << nx << " double\n";
            for (auto i: range(0,nx)) out_str << coord_sys.xcoord.map(bnd.min(0)+i*obj.get_dx(0)) << "\n";
            out_str << "Y_COORDINATES " << ny << " double\n";
            for (auto j: range(0,ny)) out_str << coord_sys.ycoord.map(bnd.min(1)+j*obj.get_dx(1)) << "\n";
            out_str << "Z_COORDINATES " << nz << " double\n";
            for (auto k: range(0,nz)) out_str << coord_sys.zcoord.map(bnd.min(2)+k*obj.get_dx(2)) << "\n";
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, coords::dense_coordinate_system dense_coord_sys_t>
        static void output_mesh_data(output_stream_t& out_str, const grid_output_t& obj, const dense_coord_sys_t& coord_sys)
        {
            out_str << "DATASET STRUCTURED_GRID\nDIMENSIONS ";
            ctrs::array<int, 3> dims(
                1+obj.get_num_blocks(0)*obj.get_num_cells(0),
                1+obj.get_num_blocks(1)*obj.get_num_cells(1),
                (obj.dim()==3)*(1)+obj.get_num_blocks(2)*obj.get_num_cells(2));
            out_str << dims[0] <<  " ";
            out_str << dims[1] <<  " ";
            out_str << dims[2] << "\n";
            out_str << "POINTS " << dims[0]*dims[1]*dims[2] << " double\n";
            
            auto r0 = range(0,obj.get_num_cells(0))*range(0, obj.get_num_blocks(0));
            auto r1 = range(0,obj.get_num_cells(1))*range(0, obj.get_num_blocks(1));
            auto r2 = range(0,obj.get_num_cells(2))*range(0, obj.get_num_blocks(2));

            int i, j, k, lbi, lbj, lbk, lb;
            int nlbi = obj.get_num_blocks(0);
            int nlbj = obj.get_num_blocks(1);
            int nlbk = obj.get_num_blocks(2);
            lbj = 0;
            lbk = 0;
            j = 0;
            k = 0;
            std::size_t ct = 0;
            ctrs::array<typename dense_coord_sys_t::coord_type, 3> x;
            for (lbk = 0; lbk < obj.get_num_blocks(2); ++lbk)
            {
                for (k = 0; k < obj.get_num_cells(2)+((lbk==(nlbk-1))?(1):(0)); ++k)
                {
                    for (lbj = 0; lbj < obj.get_num_blocks(1); ++lbj)
                    {
                        for (j = 0; j < obj.get_num_cells(1)+((lbj==(nlbj-1))?(1):(0)); ++j)
                        {
                            for (lbi = 0; lbi < obj.get_num_blocks(0); ++lbi)
                            {
                                lb = obj.collapse_block_num(lbi,lbj,lbk);
                                for (i = 0; i < obj.get_num_cells(0)+((lbi==(nlbi-1))?(1):(0)); ++i)
                                {
                                    grid::node_idx_t i_n(i, j, k, lb);
                                    x = obj.get_coords(i_n);
                                    out_str << x[0] << " " << x[1] << " " << x[2] << "\n";
                                    ct++;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, typename... arrays_t>
        static void output_grid_serial(output_stream_t& out_str, const grid_output_t& obj, arrays_t... arrays)
        {
            
            output_serial_header(out_str);
            output_mesh_data(out_str, obj, obj.coord_sys());
            
            auto test_node_center = [](auto i) -> std::size_t {return (i.centering_type()==grid::node_centered)?1:0;};
            std::size_t num_node_centered_arrays = utils::reduce_over_params(test_node_center, arrays...);
            
            auto test_cell_center = [](auto i) -> std::size_t {return (i.centering_type()==grid::cell_centered)?1:0;};
            std::size_t num_cell_centered_arrays = utils::reduce_over_params(test_cell_center, arrays...);
            
            int ct = 0;
            auto current_centering = grid::node_centered;
            auto output_array_data = [&](auto& arr) -> void
            {
                if (arr.centering_type()==current_centering)
                {
                    for (auto j: range(0,arr.get_major_dims().total_size()))
                    {
                        for (auto i: range(0,arr.get_minor_dims().total_size()))
                        {
                            std::string name = get_array_name(arr, i, j);
                            out_str << "SCALARS " << name << " double\nLOOKUP_TABLE default\n";
                            auto grid_range = obj.get_range(arr.centering_type(), grid::exclude_exchanges);
                            auto r1 = grid_range.subrange(0)*range(0, obj.get_num_blocks(0));
                            auto r2 = grid_range.subrange(1)*range(0, obj.get_num_blocks(1));
                            auto r3 = grid_range.subrange(2)*range(0, obj.get_num_blocks(2));
                            auto total_range = r1*r2*r3;
                            
                            for (auto midx: total_range)
                            {
                                std::size_t lb = obj.collapse_block_num(midx[1], midx[3], midx[5]);
                                out_str << arr.unwrap_idx(i, midx[0], midx[2], midx[4], lb, j) << "\n";
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
        
        static std::string ntab(const std::size_t& n, const std::size_t& tab_size = 4) { return std::string(tab_size*n, ' '); }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t>
        static void output_parallel_header_file(
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
                out_str << ntab(3) << utils::strformat("<DataSet index=\"{}\" file=\"{}\"/>\n", i, utils::strformat(block_proto_string, utils::zfill(i,num_zeros)));
            }
            out_str << ntab(2) << "</Block>\n";
            out_str << ntab(1) << "</vtkMultiBlockDataSet>\n";
            out_str << "</VTKFile>\n";
        }
        
        template <class output_stream_t, grid::multiblock_grid grid_output_t, typename... arrays_t>
        static void output_parralel_block_file(output_stream_t& out_str, const std::size_t& lb_loc, const grid_output_t& obj, arrays_t... arrays)
        {
            out_str << "<?xml version=\"1.0\"?>\n";
            out_str << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"" << (utils::is_big_endian()?"BigEndian":"LittleEndian") << "\" header_type=\"UInt32\">\n";
            const int n_total_i = obj.get_num_cells(0)+2*obj.get_num_exchange(0);
            const int n_total_j = obj.get_num_cells(1)+2*obj.get_num_exchange(1);
            const int n_total_k = obj.get_num_cells(2)+2*obj.get_num_exchange(2);
            const int n_cells_i = obj.get_num_cells(0);
            const int n_cells_j = obj.get_num_cells(1);
            const int n_cells_k = obj.get_num_cells(2);
            const int n_guard_i = obj.get_num_exchange(0);
            const int n_guard_j = obj.get_num_exchange(1);
            const int n_guard_k = obj.get_num_exchange(2);
            auto box = obj.get_block_box(lb_loc);
            out_str << ntab(1) << utils::strformat("<RectilinearGrid WholeExtent=\"0 {} 0 {} 0 {}\">", n_total_i, n_total_j, n_total_k) << std::endl;
            out_str << ntab(2) << "<FieldData>" << std::endl;
            out_str << ntab(3) << "<DataArray type=\"Int32\" Name=\"avtRealDims\" NumberOfTuples=\"6\" format=\"ascii\">" << std::endl;
            out_str << ntab(4) << utils::strformat("{} {} {} {} {} {}", n_guard_i, n_guard_i+n_cells_i, n_guard_j, n_guard_j+n_cells_j, n_guard_k, obj.is_3d()?(n_guard_k+n_cells_k):0) << std::endl;
            out_str << ntab(3) << "</DataArray>" << std::endl;
            out_str << ntab(3) << "<DataArray type=\"Float64\" Name=\"avtOriginalBounds\" NumberOfTuples=\"6\" format=\"ascii\">" << std::endl;
            out_str << ntab(4) << utils::strformat("{} {} {} {} {} {}", box.min(0), box.max(0), box.min(1), box.max(1), box.min(2), box.max(2)) << std::endl;
            out_str << ntab(3) << "</DataArray>" << std::endl;
            out_str << ntab(2) << "</FieldData>" << std::endl;
            out_str << ntab(2) << utils::strformat("<Piece Extent=\"0 {} 0 {} 0 {}\">", n_total_i, n_total_j, n_total_k) << std::endl;
            
            std::string vars_string = "";
            std::size_t idx = 0;
            auto add_name = [&](auto& i) -> void
            {
                for (auto i1:range(0,i.get_minor_dims().total_size()))
                {
                    for (auto i2:range(0,i.get_major_dims().total_size()))
                    {
                        if (idx>0) vars_string+=",";
                        vars_string += get_array_name(i, i1, i2);
                        ++idx;
                    }
                }
            };
            utils::foreach_param(add_name, arrays...);
            out_str << ntab(3) << "<CellData Scalars=\"" << vars_string << "\">" << std::endl;
            auto write_data = [&](const auto& arr) -> void
            {
                
                typedef typename std::remove_reference<decltype(arr)>::type::value_type data_t;
                
                if (arr.centering_type()!=grid::cell_centered) throw std::runtime_error("parallel output not currently supporting anythong other than cell data!");
                std::size_t total_temp_size = (n_cells_i+2*n_guard_i)*(n_cells_j+2*n_guard_j)*(n_cells_k+2*n_guard_k);
                unsigned int total_bytes_size = sizeof(data_t)*(unsigned int)total_temp_size;
                std::vector<data_t> compressed_data(total_temp_size);
                auto block_grid_range = range(-n_guard_i, n_cells_i + n_guard_i)*range(-n_guard_j, n_cells_j + n_guard_j)*range(-n_guard_k, n_cells_k + n_guard_k);
                for (auto i2:range(0,arr.get_major_dims().total_size()))
                {
                    for (auto i1:range(0,arr.get_minor_dims().total_size()))
                    {
                        std::size_t idx = 0;
                        for (auto ijk: block_grid_range)
                        {
                            compressed_data[idx++] = arr.unwrap_idx(i1, ijk[0], ijk[1], ijk[2], lb_loc, i2);
                        }
                        out_str << ntab(4) << "<DataArray type=\"Float64\" Name=\"" << get_array_name(arr, i1, i2) << "\" format=\"binary\">\n";
                        //data here
                        // here here he
                        spade::detail::stream_base_64(out_str, &total_bytes_size, 1);
                        spade::detail::stream_base_64(out_str, &compressed_data[0], compressed_data.size());
                        //do this bit right now
                        out_str << "\n" << ntab(4) << "</DataArray>\n";
                    }
                }
            };
            utils::foreach_param(write_data, arrays...);
            out_str << ntab(4) << "<DataArray type=\"UInt8\" Name=\"avtGhostZones\" format=\"ascii\">\n";
            auto isGhost = [&](int i, int j, int k) -> bool
            {
                return (i<n_guard_i)||(i>=n_total_i-n_guard_i)||(j<n_guard_j)||(j>=n_total_j-n_guard_j)||(k<n_guard_k)||(k>=n_total_k-n_guard_k);
            };
            ctrs::array<double,3> xyz(0, 0, 0);
            std::string csp20 = ntab(5);
            for (auto ijk: range(0,n_total_i)*range(0,n_total_j)*range(0,n_total_k))
            {
                out_str << csp20 << (isGhost(ijk[0], ijk[1], ijk[2])?16:0) << "\n";
            }
            out_str << ntab(4) << "</DataArray>\n";
            out_str << ntab(3) << "</CellData>\n";
            out_str << ntab(3) << "<Coordinates>" << std::endl;
            
            ctrs::array<int,3> ng(n_guard_i, n_guard_j, n_guard_k);
            auto box_tmp = obj.get_block_box(lb_loc);
            algs::static_for<0, obj.dim()>([&](auto i)->void
            {
                box_tmp.min(i.value) -= obj.get_dx(i.value)*ng[i.value];
                box_tmp.max(i.value) += obj.get_dx(i.value)*ng[i.value];
            });
            
            out_str << ntab(4) << spade::utils::strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", box_tmp.min(0), box_tmp.max(0)) << std::endl;
            int i,j,k;
            j = 0;
            k = 0;
            for (i = -n_guard_i; i <=n_cells_i+n_guard_i; i++)
            {
                grid::node_idx_t i_n(i, j, k, lb_loc);
                xyz = obj.get_coords(i_n);
                out_str << csp20 << xyz[0] << "\n";
            }
            out_str << ntab(4) << "</DataArray>" << std::endl;
            out_str << ntab(4) << spade::utils::strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", box_tmp.min(1), box_tmp.max(1)) << std::endl;
            
            i = 0;
            k = 0;
            for (j = -n_guard_j; j <=n_cells_j+n_guard_j; j++)
            {
                grid::node_idx_t i_n(i, j, k, lb_loc);
                xyz = obj.get_coords(i_n);
                out_str << csp20 << xyz[1] << "\n";
            }
            out_str << ntab(4) << "</DataArray>" << std::endl;
            out_str << ntab(4) << spade::utils::strformat("<DataArray type=\"Float64\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">", box_tmp.min(2), box_tmp.max(2)) << std::endl;
            
            i = 0;
            j = 0;
            for (k = -n_guard_k; k <=n_cells_k+n_guard_k; k++)
            {
                grid::node_idx_t i_n(i, j, k, lb_loc);
                xyz = obj.get_coords(i_n);
                out_str << csp20 << xyz[2] << "\n";
            }
            out_str << ntab(4) << "</DataArray>" << std::endl;
            out_str << ntab(3) << "</Coordinates>" << std::endl;
            out_str << ntab(2) << "</Piece>" << std::endl;
            out_str << ntab(1) << "</RectilinearGrid>" << std::endl;
            out_str << "</VTKFile>" << std::endl;
        }
        
        template <grid::multiblock_grid grid_output_t, typename... arrays_t>
        static std::string output_grid_parallel(const std::string& out_dir, const std::string& out_name_no_extension, const grid_output_t& obj, arrays_t... arrays)
        {
            std::filesystem::path out_path(out_dir);
            if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
            std::filesystem::path header_filename(out_dir);
            header_filename /= (out_name_no_extension + ".vtm");
            std::size_t numzeros = 7;
            std::string blocks_dir_name = ("data_" + out_name_no_extension);
            std::string blocks_file_name_template = "block{}.vtr";
            if (obj.group().isroot())
            {
                std::filesystem::path blocks_directory(out_dir);
                blocks_directory /= blocks_dir_name;
                if (!std::filesystem::is_directory(blocks_directory)) std::filesystem::create_directory(blocks_directory);
                
                std::filesystem::path block_file_rel_template(blocks_dir_name);
                block_file_rel_template /= blocks_file_name_template;
                std::string block_template_str = block_file_rel_template;
                std::ofstream header_file_strm(header_filename);
                
                output_parallel_header_file(header_file_strm, block_template_str, numzeros, obj);
            }
            obj.group().sync();
            for (auto i: range(0, obj.get_num_local_blocks()))
            {
                std::filesystem::path block_abs_path(out_dir);
                block_abs_path /= blocks_dir_name;
                block_abs_path /= blocks_file_name_template;
                std::string block_abs_path_str = utils::strformat(block_abs_path, utils::zfill(obj.get_partition().get_global_block(i), numzeros));
                std::ofstream block_file_strm(block_abs_path_str);
                output_parralel_block_file(block_file_strm, i, obj, arrays...);
            }
            obj.group().sync();
            return header_filename;
        }
    }
    
    template <grid::multiblock_grid obj_t, typename... arrays_t>
    static std::string output_vtk(const std::string& out_dir, const std::string& out_name_no_extension, const obj_t& obj, arrays_t... arrays)
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
    
    template <grid::multiblock_array arr_t, typename... arrays_t>
    static std::string output_vtk(const std::string& out_dir, const std::string& out_name_no_extension, const arr_t& arr, arrays_t... arrays)
    {
        return output_vtk(out_dir, out_name_no_extension, arr.get_grid(), arr, arrays...);
    }
    
    template <block_config::block_configuration blocks_t, parallel::parallel_group group_t>
    static std::string output_vtk(const std::string& out_dir, const std::string& out_name_no_extension, const blocks_t& blocks, const group_t& group)
    {
        std::filesystem::path out_path(out_dir);
        if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
        out_path /= (out_name_no_extension + ".vtk");
        std::string total_filename = out_path;
        if (group.isroot())
        {
            std::ofstream out_str(total_filename);
            const std::size_t points_per_block = static_math::pow<2,blocks.dim()>::value;
            auto is_terminal_block = [](const auto& node) ->bool {return node.terminal();};
            const std::size_t num_term_blocks = blocks.get_num_blocks(is_terminal_block);
            const std::size_t num_points = points_per_block*num_term_blocks;
            out_str << "# vtk DataFile Version 2.0\nspade output\nASCII\n";
            out_str << "DATASET UNSTRUCTURED_GRID\nPOINTS " << num_points << " double\n";
            ctrs::array<typename blocks_t::coord_val_type,3> x = 0.0;
            for (auto i: range(0, blocks.get_num_blocks()))
            {
                if (is_terminal_block(*(blocks.all_nodes[i])))
                {
                    const auto bbox = blocks.get_block_box(i);
                    for (auto min_max: range(0,2)*range(0,2)*range(0,blocks.dim()==3?2:1))
                    {
                        for (auto j: range(0,blocks.dim()))
                        {
                            x[j] = bbox(j, min_max[j]);
                        }
                        out_str << x[0] << " " << x[1] << " " << x[2] << "\n";
                    }
                }
            }
            out_str << "CELLS " << num_term_blocks << " " << num_term_blocks*(1+points_per_block) << "\n";
            std::size_t ct = 0;
            for (auto i: range(0, num_term_blocks))
            {
                out_str << points_per_block;
                for (auto j: range(0,points_per_block))
                {
                    out_str << " " << ct++;
                }
                out_str << "\n";
            }
            out_str << "CELL_TYPES " << num_term_blocks << "\n";
            int cell_type = (blocks.dim()==3)?11:8;
            for (auto i: range(0,num_term_blocks)) out_str << cell_type << "\n";
        }
        return total_filename;
    }
    
    template <typename idx_t, grid::multiblock_grid grid_t>
    static std::string output_vtk(const std::string& out_dir, const std::string& out_name_no_extension, const std::vector<idx_t>& indices, const grid_t& grid)
    {
        const auto& group = grid.group();
        std::filesystem::path out_path(out_dir);
        if (!std::filesystem::is_directory(out_path)) std::filesystem::create_directory(out_path);
        out_path /= (out_name_no_extension + ".vtk");
        std::string total_filename = out_path;
        using coord_type = typename grid_t::coord_type;
        std::vector<coord_type> points;
        points.resize(3*indices.size());
        std::size_t id = 0;
        for (const auto& idx: indices)
        {
            const auto x = grid.get_coords(idx);
            points[id++] = x[0];
            points[id++] = x[1];
            points[id++] = x[2];
        }
        const std::size_t total_points = group.sum(indices.size());
        std::vector<coord_type> all_data;
        group.append_root(all_data, points);
        if (group.isroot())
        {
            std::ofstream out_str(total_filename);
            out_str << "# vtk DataFile Version 2.0\nspade output\nASCII\nDATASET POLYDATA\nPOINTS ";
            out_str << total_points << " double\n";
            std::size_t ct = 0;
            for (auto ll: range(0, total_points))
            {
                out_str << all_data[ct++] << " " << all_data[ct++] << " " << all_data[ct++] << "\n";
            }
            out_str << "POINT_DATA " << total_points << "\nSCALARS Data double\nLOOKUP_TABLE default\n";
            for (auto ll: range(0, total_points)) out_str << "1\n";
        }
        return total_filename;
    }
    
    namespace detail
    {
        template <grid::multiblock_array array_t> void get_array_par_buf(parallel::par_buf_t& buf, const array_t& array)
        {
            buf.clear();
            const auto& grid = array.get_grid();
            auto block_range = grid.get_range(array.centering_type(), grid::include_exchanges);
            std::size_t block_elems = array.get_minor_dims().total_size()*block_range.size() / grid.get_num_local_blocks();
            typedef typename array_t::value_type data_t;
            std::size_t block_size_bytes = block_elems*sizeof(data_t);
            ctrs::array<int, 3> nexch (grid.get_num_exchange(0), grid.get_num_exchange(1), grid.get_num_exchange(2));
            for (auto maj: range(0, array.get_major_dims().total_size()))
            {
                for (auto lb: range(0, grid.get_num_local_blocks()))
                {
                    std::size_t lb_glob = array.get_grid().get_partition().get_global_block(lb);
                    void* ptr = (void*)(&array.unwrap_idx(0,-nexch[0], -nexch[1], -nexch[2], lb, maj));
                    std::size_t offset_bytes = lb_glob*block_size_bytes;
                    buf.add(ptr, block_size_bytes, offset_bytes);
                }
            }
        }
    }
    
    template <grid::multiblock_array array_t> void binary_write(const std::string& filename, array_t& array)
    {
        const auto& group = array.get_grid().group();
        parallel::par_buf_t buf;
        detail::get_array_par_buf(buf, array);
        parallel::mpi_file_t mf(group, filename);
        mf.write_buf(buf);
    }
    
    template <grid::multiblock_array array_t> void binary_read(const std::string& filename, array_t& array)
    {
        const auto& group = array.get_grid().group();
        parallel::par_buf_t buf;
        detail::get_array_par_buf(buf, array);
        parallel::mpi_file_t mf(group, filename);
        mf.read_buf(buf);
    }
}