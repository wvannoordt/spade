#pragma once

#include <fstream>
#include <filesystem>

#include "core/utils.h"
#include "core/except.h"
#include "core/base_64.h"
#include "grid/grid.h"
#include "dispatch/execute.h"
#include "dispatch/support_of.h"

#include "ibm/boundary_info.h"

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
        
        template <typename data_t, typename idx_t>
        requires (ctrs::basic_array<data_t>)
        _sp_hybrid static auto ith_elem(const data_t& d, const idx_t& i) { return d[i]; }
        
        template <typename data_t, typename idx_t>
        _sp_hybrid static auto ith_elem(const data_t& d, const idx_t&) { return d; }
        
        template <typename arr_t, typename vec_t, typename var_idx_t>
        static void copy_block_variable_data(
            const arr_t& arr,
            vec_t& data,
            const var_idx_t vidx,
            const std::size_t& lb)
        {
            constexpr bool is_gpu_vec = requires { data.devc_data; };
            
            auto range = dispatch::support_of(arr, grid::exclude_exchanges);
            range.lower.lb() = lb;
            range.upper.lb() = lb + 1;
            
            std::size_t ni = arr.get_grid().get_num_cells(0);
            std::size_t nj = arr.get_grid().get_num_cells(1);
            std::size_t nk = arr.get_grid().get_num_cells(2);
            
            using data_t = vec_t::value_type;
            data_t* ptr = &data[0];
            if constexpr (is_gpu_vec) ptr = &(data.devc_data[0]);
            auto img = arr.image();
            using index_type = decltype(range)::index_type;
            auto load = _sp_lambda (const index_type& i) mutable
            {
                std::size_t offst = i.i() + ni*i.j() + ni*nj*i.k();
                const auto data   = img.get_elem(i);
                ptr[offst] = ith_elem(data, vidx);
            };
            
            dispatch::execute(range, load);
            
            if constexpr (is_gpu_vec) data.itransfer();
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
            bf << "\" GhostLevel=\"" << utils::min(arr.get_num_exchange(0), arr.get_num_exchange(1), arr.get_num_exchange(2), 0) << "\">\n";
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
            const auto& geom  = grid.geometry(partition::local);
            std::filesystem::path out_path(out_dir);
            using coord_float_t = utils::remove_all<decltype(grid)>::type::coord_type;
            using data_float_t  = typename arr_t::value_type;
            using alias_type    = typename arr_t::alias_type;
            const std::string data_direc = out_path / get_data_basename(basename);
            if (!std::filesystem::is_directory(data_direc) && group.isroot()) std::filesystem::create_directory(data_direc);
            group.sync();
            const auto bbx = grid.get_bounds();
            bound_box_t<std::size_t, 3> vtk_extent;
            vtk_extent.min(0) = 0;
            vtk_extent.min(1) = 0;
            vtk_extent.min(2) = 0;
            vtk_extent.max(0) = grid.get_num_cells(0);// - 1;
            vtk_extent.max(1) = grid.get_num_cells(1);// - 1;
            vtk_extent.max(2) = grid.get_num_cells(2);// - 1;
            std::vector<coord_float_t> coord_raw;
            
            using device_t = typename arr_t::device_type;
            constexpr bool is_gpu_device = device::is_gpu<device_t>;
            using vec_t = std::conditional<is_gpu_device, device::shared_vector<data_float_t, device::pinned_allocator_t<data_float_t>, device::device_allocator_t<data_float_t>>, std::vector<data_float_t>>::type;
            vec_t data_raw;
            
            coord_raw.reserve(3*(grid.get_num_cells(0)+1)*(grid.get_num_cells(1)+1)*(grid.get_num_cells(2)+1));
            data_raw.resize(grid.get_num_cells(0)*grid.get_num_cells(1)*grid.get_num_cells(2));
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
                    copy_block_variable_data(arr, data_raw, ct, lb);
                    
                    bf << "<DataArray type=\"" << get_fundamental_str(data_float_t()) << "\" Name=\"" << name << "\" format=\"" << format_str << "\">\n";
                    const auto& rdat = [&]()
                    {
                        if constexpr (is_gpu_device)
                        {
                            return data_raw.host_data;
                        }
                        else
                        {
                            return data_raw;
                        }
                    }();
                    spade::detail::stream_base_64(bf, rdat);
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
                            auto pt = geom.get_coords(idx);
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
        group.sync();
        detail::output_block_files(arr, out_dir, basename);
        group.sync();
    }
    
    template <ctrs::basic_array ctr_t>
    requires(std::same_as<typename ctr_t::value_type, coords::point_t<typename ctr_t::value_type::value_type>>)
    static void output_vtk(const std::string& fname, const ctr_t& data, const bool fill_line = false)
    {
        std::ofstream mf(fname);
        mf << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << data.size() << " double\n";
        for (const auto& x: data) mf << x[0] << " " << x[1] << " " << x[2] << "\n";
        if (fill_line && (data.size() > 1))
        {
            mf << "POLYGONS " << (data.size() - 1) << " " << 3*(data.size() - 1) << "\n";
            for (std::size_t i = 0; i < data.size()-1; ++i)
            {
                mf << "2 " << i << " " << (i+1) << "\n";
            }
        }
        mf << "POINT_DATA " << data.size() << "\nSCALARS Data double\nLOOKUP_TABLE default\n";
        std::size_t ct = 0;
        for (const auto& x: data) mf << ct++ << "\n";
    }
    
    namespace detail
    {
        template <std::integral int_t>
        inline std::string get_ascii_name(const int_t&) { return "int"; }
        
        template <std::floating_point dbl_t>
        inline std::string get_ascii_name(const dbl_t&) { return "double"; }
        
        inline void cell_data_out_ascii(std::ostream& mf, std::size_t required_size) {}
        
        template <typename data_t, typename... others_t>
        inline void cell_data_out_ascii(std::ostream& mf, std::size_t required_size, const std::string& scal_name, const std::vector<data_t>& data, const others_t&... others)
        {
            if (data.size() != required_size)
            {
                throw except::sp_exception("attempted to write to VTK file with incorrectly-sized scalar \"" + scal_name + "\"!");
            }
            mf << "SCALARS " << scal_name << " " << get_ascii_name(data_t()) << "\n";
            mf << "LOOKUP_TABLE default\n";
            for (const auto& d: data) mf << d << "\n";
            cell_data_out_ascii(mf, required_size, others...);
        }
    }
    
    template <typename bbx_ctr_t, typename... others_t>
    requires (requires { bbx_ctr_t()[0].min(0); bbx_ctr_t()[0].max(0); bbx_ctr_t()[0].size(0); }) // good lord
    static void output_vtk(const std::string& fname, const bbx_ctr_t& blist, const others_t&... others)
    {
        constexpr static int bdim = bbx_ctr_t::value_type::size();
        std::ofstream mf(fname);
        static_assert((bdim == 3) || (bdim == 2), "io only supports 2d/3d bounding box");
        std::size_t npts = 8*blist.size();
        if constexpr (bdim==2) npts = 4*blist.size();
        mf << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << npts << " double\n";
        for (std::size_t i = 0; i < blist.size(); ++i)
        {
            const auto& bnd = blist[i];
            if constexpr (bdim == 3)
            {
                mf << bnd.min(0) << " " << bnd.min(1) << " " << bnd.min(2) << "\n";
                mf << bnd.max(0) << " " << bnd.min(1) << " " << bnd.min(2) << "\n";
                mf << bnd.min(0) << " " << bnd.max(1) << " " << bnd.min(2) << "\n";
                mf << bnd.max(0) << " " << bnd.max(1) << " " << bnd.min(2) << "\n";
                mf << bnd.min(0) << " " << bnd.min(1) << " " << bnd.max(2) << "\n";
                mf << bnd.max(0) << " " << bnd.min(1) << " " << bnd.max(2) << "\n";
                mf << bnd.min(0) << " " << bnd.max(1) << " " << bnd.max(2) << "\n";
                mf << bnd.max(0) << " " << bnd.max(1) << " " << bnd.max(2) << "\n";
            }
            else
            {
                mf << bnd.min(0) << " " << bnd.min(1) << "\n";
                mf << bnd.max(0) << " " << bnd.min(1) << "\n";
                mf << bnd.min(0) << " " << bnd.max(1) << "\n";
                mf << bnd.max(0) << " " << bnd.max(1) << "\n";
            }
        }

        if constexpr (bdim == 3) mf << "CELLS " << blist.size() << " " << 9*blist.size() << "\n";
        if constexpr (bdim == 2) mf << "CELLS " << blist.size() << " " << 5*blist.size() << "\n";
        std::size_t ipt = 0;
        for (std::size_t i = 0; i < blist.size(); ++i)
        {
            mf << "8";
            for (int pp = 0; pp < ((bdim==3)?8:4); ++pp)
            {
                mf << " " << ipt++;
            }
            mf << "\n";
        }
        mf << "CELL_TYPES " << blist.size() << "\n";
        for (std::size_t pp = 0; pp < blist.size(); ++pp)
        {
            mf << ((bdim==3)?11:8) << "\n";
        }
        
        if constexpr(sizeof...(others_t) > 0)
        {
            mf << "CELL_DATA " << blist.size() << "\n";
            detail::cell_data_out_ascii(mf, blist.size(), others...);
        }
    }
    
    template <grid::multiblock_grid grid_t>
    static void output_vtk(const std::string& fname, const grid_t& grid)
    {
        const auto& group = grid.group();
        using bbx_type = decltype(grid.get_bounding_box(0));
        std::vector<bbx_type> boxes;
        std::vector<int> ranks;
        std::vector<int> block_local;
        std::vector<int> block_global;
        
        const auto& pp = grid.get_partition();
        
        for (std::size_t lb = 0; lb < grid.get_num_global_blocks(); ++lb)
        {
            auto lbt = utils::tag[partition::global](lb);
            ranks.push_back(pp.get_rank(lbt));
            block_local.push_back(pp.get_any_local(lbt)); //TODO
            block_global.push_back(lb);
            boxes.push_back(grid.get_bounding_box(lbt));
        }
        
        if (group.isroot()) output_vtk(fname, boxes, "rank", ranks, "block_id_local", block_local, "block_id_global", block_global);
    }
    
    template <
        const int grid_dim,
        const int nlayers,
        typename float_t,
        template <typename> typename container_t,
        grid::multiblock_grid grid_t>
    inline void output_vtk(
        const std::string& direc,
        const std::string& title,
        const ibm::boundary_info_t<grid_dim, nlayers, float_t, container_t>& info,
        const grid_t& grid)
    {
        
        const auto make_filename = [&](const std::string& thing)
        {
            std::filesystem::path out_path = direc;
            out_path /= title;
            std::string output = out_path;
            output += ".";
            output += thing;
            output += ".vtk";
            return output;
        };
        
        std::ofstream mf(make_filename("ghosts"));
        ctrs::array<std::size_t, grid_t::dim()> n_elems = 0;
        std::size_t total_points = 0;
        for (int d = 0; d < grid_t::dim(); ++d)
        {
            n_elems[d]    = info.aligned[d].indices.size();
            total_points += n_elems[d]*nlayers;
        }
        
        mf << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << total_points << " double\n";
        for (int d = 0; d < grid_t::dim(); ++d)
        {
            for (std::size_t ii = 0; ii < n_elems[d]; ++ii)
            {
                for (int lyr = 0; lyr < nlayers; ++lyr)
                {
                    const auto icell = info.aligned[d].indices[ii][lyr];
                    const auto xx = grid.get_coords(icell);
                    mf << xx[0] << " " << xx[1] << " " << xx[2] << "\n";
                }
            }
        }
        
        mf << "POINT_DATA " << total_points << "\n";
        
        const auto scalar_out = [&](auto& flh, const std::string& scname, const auto& fnc)
        {
            flh << "SCALARS " << scname << " double\nLOOKUP_TABLE default\n";
            for (int d = 0; d < grid_t::dim(); ++d)
            {
                for (std::size_t ii = 0; ii < n_elems[d]; ++ii)
                {
                    for (int lyr = 0; lyr < nlayers; ++lyr)
                    {
                        flh << fnc(d, ii, lyr) << "\n";
                    }
                }
            }
        };
        
        const auto vector_out = [&](auto& flh, const std::string& scname, const auto& fnc)
        {
            flh << "VECTORS " << scname << " double\n";
            for (int d = 0; d < grid_t::dim(); ++d)
            {
                for (std::size_t ii = 0; ii < n_elems[d]; ++ii)
                {
                    for (int lyr = 0; lyr < nlayers; ++lyr)
                    {
                        auto vl = fnc(d, ii, lyr);
                        flh << vl[0] << " " << vl[1] << " " << vl[2] << "\n";
                    }
                }
            }   
        };
        
        scalar_out(mf, "layer",       [&](const int dir, const int id, const int layer){ return layer; });
        scalar_out(mf, "fillable",    [&](const int dir, const int id, const int layer){ return info.aligned[dir].can_fill[id][layer]; });
        scalar_out(mf, "direction",   [&](const int dir, const int id, const int layer){ return dir; });
        scalar_out(mf, "id",          [&](const int dir, const int id, const int layer){ return id; });
        scalar_out(mf, "i",           [&](const int dir, const int id, const int layer){ return info.aligned[dir].indices[id][layer].i();  });
        scalar_out(mf, "j",           [&](const int dir, const int id, const int layer){ return info.aligned[dir].indices[id][layer].j();  });
        scalar_out(mf, "k",           [&](const int dir, const int id, const int layer){ return info.aligned[dir].indices[id][layer].k();  });
        scalar_out(mf, "lb",          [&](const int dir, const int id, const int layer){ return info.aligned[dir].indices[id][layer].lb(); });
        vector_out(mf, "boundary_pt", [&](const int dir, const int id, const int layer){ return info.aligned[dir].boundary_points[id]       - grid.get_coords(info.aligned[dir].indices[id][layer]); });
        vector_out(mf, "closest_pt",  [&](const int dir, const int id, const int layer){ return info.aligned[dir].closest_points[id][layer] - grid.get_coords(info.aligned[dir].indices[id][layer]); });
        
        
        std::ofstream mf2(make_filename("bndy_pts"));
        mf2 << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << total_points << " double\n";
        for (int d = 0; d < grid_t::dim(); ++d)
        {
            for (std::size_t ii = 0; ii < n_elems[d]; ++ii)
            {
                for (int lyr = 0; lyr < nlayers; ++lyr)
                {
                    const auto xx = info.aligned[d].closest_points[ii][lyr];
                    mf2 << xx[0] << " " << xx[1] << " " << xx[2] << "\n";
                }
            }
        }
        mf2 << "POINT_DATA " << total_points << "\n";
        scalar_out(mf2, "layer",       [&](const int dir, const int id, const int layer){ return layer; });
        scalar_out(mf2, "direction",   [&](const int dir, const int id, const int layer){ return dir;   });
    }
    
}
