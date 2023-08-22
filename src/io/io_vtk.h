#pragma once

#include <fstream>
#include <filesystem>

#include "core/utils.h"
#include "core/except.h"
#include "core/base_64.h"

namespace spade::io
{
    // Resources:
    // https://examples.vtk.org/site/VTKFileFormats/
    // https://discourse.vtk.org/t/creating-a-parallel-vtk-structured-grid/9009
    namespace detail
    {
        static std::string endian_str()
        {
            return utils::is_big_endian()?"BigEndian":"LittleEndian";
        }
        
        template <typename T> concept named_array = requires(T t) { t.name(T::size() - 1); };
        template <typename T> concept unnamed_array = ctrs::basic_array<T> && !named_array<T>;
        
        template <named_array alias_t>
        std::vector<std::string> get_var_names(const alias_t&)
        {
            std::vector<std::string> pp(alias_t::size());
            for (int i = 0; i < pp.size(); ++i) pp[i] = alias_t::name(i);
            return pp;
        }
        
        template <unnamed_array alias_t>
        std::vector<std::string> get_var_names(const alias_t&)
        {
            std::vector<std::string> pp(alias_t::size(), "data");
            for (int i = 0; i < pp.size(); ++i) pp[i] += utils::zfill(i, 3);
            return pp;
        }
        
        template <typename alias_t>
        requires (!named_array<alias_t> && !unnamed_array<alias_t>)
        std::vector<std::string> get_var_names(const alias_t&)
        {
            std::vector<std::string> pp{"data"};
            return pp;
        }
        
        template <typename fundamental_t>
        static std::string get_fundamental_str(const fundamental_t&)
        {
            if constexpr (std::same_as<fundamental_t, int>)    return "Int32";
            if constexpr (std::same_as<fundamental_t, float>)  return "Float32";
            if constexpr (std::same_as<fundamental_t, double>) return "Float64";
            throw except::sp_exception("unsupported fundamental type in vtk output");
            return "[ERROR]";
        }
        
        static std::string get_data_basename(const std::string& base_name)
        {
            return "data_" + base_name;
        }
        
        static std::string get_block_filename(const std::string& base_name, const std::size_t& lb_glob)
        {
            std::filesystem::path of(get_data_basename(base_name));
            of /= (std::string("b") + utils::zfill(lb_glob, 9) + ".vts");
            return std::string(of);
        }
        
        static std::string get_center_str(const grid::array_centering& i)
        {
            switch(i)
            {
                case grid::cell_centered:
                {
                    return "CellData";
                }
                default:
                {
                    throw except::sp_exception("attempt to output unsupported variable centering type");
                }
            }
        }
        
        template <typename arr_t>
        static void output_base_file(const arr_t& arr, const std::string& base_file, const std::string& base_name)
        {
            const auto& grid    = arr.get_grid();
            using coord_float_t = utils::remove_all<decltype(grid)>::type::coord_type;
            const auto bbx      = grid.get_bounds();
            bound_box_t<std::size_t, 3> vtk_extent;
            vtk_extent.min(0) = 0;
            vtk_extent.min(1) = 0;
            vtk_extent.min(2) = 0;
            vtk_extent.max(0) = grid.get_num_cells(0);// - 1;
            vtk_extent.max(1) = grid.get_num_cells(1);// - 1;
            vtk_extent.max(2) = grid.get_num_cells(2) * grid.get_num_global_blocks();// - 1;
            const auto varnames = get_var_names(typename arr_t::alias_type());
            std::string data_str = "P" + get_center_str(arr.centering_type());
            const std::string fundamental = get_fundamental_str(typename arr_t::value_type());
            std::ofstream bf(base_file);
            bf << "<?xml version=\"1.0\"?>\n";
            bf << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"" << endian_str() << "\">\n";
            bf << "<PStructuredGrid WholeExtent=\"";
            bf << utils::join(vtk_extent.data(), " ");
            bf << "\" GhostLevel=\"" << utils::min(grid.get_num_exchange(0), grid.get_num_exchange(1), grid.get_num_exchange(2), 0) << "\">\n";
            bf << "<" << data_str << " Scalars = \"" << utils::join(varnames, ",") << "\">\n";
            for (const auto& name: varnames)
            {
                bf << "<PDataArray type=\"" << fundamental << "\" Name=\"" << name << "\"/>\n";
            }
            bf << "</" << data_str << ">\n";
            bf << "<PPoints>\n";
            bf << "<PDataArray type=\"" << get_fundamental_str(coord_float_t()) << "\" ";
            bf << "NumberOfComponents=\"3\"/>\n";
            bf << "</PPoints>\n";
            for (std::size_t lb = 0; lb < grid.get_num_global_blocks(); ++lb)
            {
                const std::string lbfn = get_block_filename(base_name, lb);
                bound_box_t<std::size_t, 3> blk_extent = vtk_extent;
                vtk_extent.min(2) = lb     * grid.get_num_cells(2);
                vtk_extent.max(2) = (lb+1) * grid.get_num_cells(2);// - 1;
                bf << "<Piece Extent=\"" << utils::join(vtk_extent.data(), " ") << "\" Source=\"" << lbfn << "\"/>\n";
            }
            bf << "</PStructuredGrid>\n";
            bf << "</VTKFile>\n";
        }
        
        template <typename arr_t>
        static void output_block_files(const arr_t& arr, const std::string& out_dir, const std::string& basename)
        {
            const auto& grid  = arr.get_grid();
            const auto& group = grid.group();
            std::filesystem::path out_path(out_dir);
            using coord_float_t = utils::remove_all<decltype(grid)>::type::coord_type;
            using data_float_t  = typename arr_t::value_type;
            using alias_type    = typename arr_t::alias_type;
            const std::string data_direc = out_path / get_data_basename(basename);
            if (!std::filesystem::is_directory(data_direc) && group.isroot()) std::filesystem::create_directory(data_direc);
            const auto bbx = grid.get_bounds();
            bound_box_t<std::size_t, 3> vtk_extent;
            vtk_extent.min(0) = 0;
            vtk_extent.min(1) = 0;
            vtk_extent.min(2) = 0;
            vtk_extent.max(0) = grid.get_num_cells(0);// - 1;
            vtk_extent.max(1) = grid.get_num_cells(1);// - 1;
            vtk_extent.max(2) = grid.get_num_cells(2);// - 1;
            std::vector<coord_float_t> coord_raw;
            std::vector<data_float_t> data_raw;
            coord_raw.reserve(3*(grid.get_num_cells(0)+1)*(grid.get_num_cells(1)+1)*(grid.get_num_cells(2)+1));
            data_raw.reserve(grid.get_num_cells(0)*grid.get_num_cells(1)*grid.get_num_cells(2));
            const std::string format_str = "binary";
            const auto names = get_var_names(typename arr_t::alias_type());
            for (std::size_t lb = 0; lb < grid.get_num_local_blocks(); ++lb)
            {
                auto lb_glob = grid.get_partition().to_global(utils::tag[partition::local](lb));
                bound_box_t<std::size_t, 3> blk_extent = vtk_extent;
                blk_extent.min(2) = 0*lb_glob.value     * grid.get_num_cells(2);
                blk_extent.max(2) = (0*lb_glob.value+1) * grid.get_num_cells(2);// - 1;
                const std::string lbfn = std::filesystem::path(out_dir) / get_block_filename(basename, lb_glob.value);
                std::ofstream bf(lbfn);
                bf << "<?xml version=\"1.0\"?>\n";
                bf << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"" << endian_str() << "\">\n";
                bf << "<StructuredGrid WholeExtent=\"" << utils::join(vtk_extent.data(), " ") << "\">\n";
                bf << "<Piece Extent=\"" << utils::join(blk_extent.data(), " ") << "\">\n";
                bf << "<"  << get_center_str(arr.centering_type()) << " Scalars=\"" << utils::join(names, ",") << "\">\n";
                int ct = 0;
                for (const auto& name: names)
                {
                    data_raw.clear();
                    
                    const auto img = arr.image();
                    grid::cell_idx_t cid(0,0,0,lb);
                    for (cid.k() = 0; cid.k() < grid.get_num_cells(2); ++cid.k())
                    {
                        for (cid.j() = 0; cid.j() < grid.get_num_cells(1); ++cid.j())
                        {
                            for (cid.i() = 0; cid.i() < grid.get_num_cells(0); ++cid.i())
                            {
                                alias_type elem = img.get_elem(cid);
                                auto data = [&]()
                                {
                                    if constexpr (ctrs::basic_array<alias_type>) return elem[ct];
                                    else return elem;
                                }();
                                data_raw.push_back(data);
                            }
                        }
                    }
                    bf << "<DataArray type=\"" << get_fundamental_str(data_float_t()) << "\" Name=\"" << name << "\" format=\"" << format_str << "\">\n";
                    spade::detail::stream_base_64(bf, data_raw);
                    bf << "\n</DataArray>\n";
                    ++ct;
                }
                bf << "</" << get_center_str(arr.centering_type()) << ">\n";
                bf << "<Points>\n";
                bf << "<DataArray type=\"" << get_fundamental_str(coord_float_t()) << "\" NumberOfComponents=\"3\" format=\"" << format_str << "\">\n";
                grid::node_idx_t idx(0,0,0,lb);
                coord_raw.clear();
                for (idx.k() = 0; idx.k() <= grid.get_num_cells(2); ++idx.k())
                {
                    for (idx.j() = 0; idx.j() <= grid.get_num_cells(1); ++idx.j())
                    {
                        for (idx.i() = 0; idx.i() <= grid.get_num_cells(0); ++idx.i())
                        {
                            auto pt = grid.get_coords(idx);
                            if constexpr (grid.dim() == 2) pt[2] = 0.0;
                            coord_raw.push_back(pt[0]);
                            coord_raw.push_back(pt[1]);
                            coord_raw.push_back(pt[2]);
                        }
                    }
                }
                spade::detail::stream_base_64(bf, coord_raw);
                bf << "\n</DataArray>\n";
                bf << "</Points>\n";
                bf << "</Piece>\n";
                bf << "</StructuredGrid>\n";
                bf << "</VTKFile>\n";
            }
        }
    }
    
    template <typename array_t>
    static void output_vtk(const std::string& out_dir, const std::string& basename, const array_t& arr)
    {
        const auto& grid  = arr.get_grid();
        const auto& group = grid.group();
        std::filesystem::path out_path(out_dir);
        if (!std::filesystem::is_directory(out_path) && group.isroot()) std::filesystem::create_directory(out_path);
        const std::string base_file  = out_path / (basename + ".pvts");
        if (group.isroot()) detail::output_base_file(arr, base_file, basename);
        detail::output_block_files(arr, out_dir, basename);
        group.sync();
    }
}
