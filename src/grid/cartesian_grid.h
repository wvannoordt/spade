#pragma once

#include "core/vec_image.h"

namespace spade::grid
{
    template <class T> concept multiblock_grid = requires(T t, const cell_idx_t& i_c, const face_idx_t& i_f, const node_idx_t& i_n)
    {
        // todo: write this
        { t.get_coords(i_c) } -> ctrs::basic_array;
        { t.get_coords(i_f) } -> ctrs::basic_array;
        { t.get_coords(i_n) } -> ctrs::basic_array;
        { t.get_comp_coords(i_c) } -> ctrs::basic_array;
        { t.get_comp_coords(i_f) } -> ctrs::basic_array;
        { t.get_comp_coords(i_n) } -> ctrs::basic_array;
    };
    
    template <class T, const array_centering ct> concept has_centering_type = (T::centering_type() == ct);
    
    template <class T> concept multiblock_array = requires(T t, int a, int i, int j, int k, int lb, int b)
    {
        t;
        t.var_map();
        t.block_map();
    };
    
    template <typename T0, typename T1> concept elementwise_compatible = multiblock_array<T0> && multiblock_array<T1> &&
    requires(const T0& t0, const T1& t1)
    {
        //TODO: write this
        t0;
    };
    
    enum exchange_inclusion_e
    {
        exclude_exchanges=0,
        include_exchanges=1
    };
    
    struct grid_dependent_t
    {
        virtual void on_blocks_update() = 0;
    };
    
    template
    <
        ctrs::basic_array array_descriptor_t,
        coords::coordinate_system coord_t,
        typename block_arrangement_t,
        parallel::parallel_group par_group_t
    >
    class cartesian_grid_t
    {
        public:
            
            using dtype              = coord_t::coord_type;
            using coord_type         = coord_t::coord_type;
            using coord_sys_type     = coord_t;
            using coord_point_type   = coords::point_t<coord_type>;
            using group_type         = par_group_t;
            using partition_type     = partition::block_partition_t;
            using dependent_type     = grid_dependent_t;
            using blocks_type        = block_arrangement_t;
            using array_desig_type   = array_descriptor_t;
            
            using geometry_type      = grid_geometry_t<coord_t, array_descriptor_t::size(), device::shared_vector>;
            using geomety_image_type = grid_geometry_t<coord_t, array_descriptor_t::size(), utils::const_vec_image_t>;
        
        private:
            partition_type grid_partition;
            geometry_type local_geometry;
            block_arrangement_t block_arrangement;
            const par_group_t& grid_group;
            std::vector<dependent_type*> dependents;

        public:
            geometry_type global_geometry;
            
            constexpr static std::size_t grid_dim = array_descriptor_t::size();
            
            cartesian_grid_t(
                const array_descriptor_t& cells_in_block_in,
                const block_arrangement_t& block_arrangement_in,
                const coord_t& coord_system_in,
                const par_group_t& group_in
                )
            : block_arrangement{block_arrangement_in},
              grid_group{group_in}
            {
                local_geometry  = geometry_type(coord_system_in, cells_in_block_in);
                global_geometry = geometry_type(coord_system_in, cells_in_block_in);
                compute_geometry();
            }
            
            //deleted until I can figure out this dependency thing
            cartesian_grid_t(const cartesian_grid_t&) = delete;
            
            cartesian_grid_t(){}            
            
            void compute_geometry()
            {
                grid_partition = partition_type(block_arrangement, grid_group);
                global_geometry.bounding_boxes.clear();
                global_geometry.domain_boundary.clear();
                local_geometry.bounding_boxes.clear();
                local_geometry.domain_boundary.clear();
                
                for (auto& v: global_geometry.boundary_blocks) v.clear();
                for (auto& v: local_geometry.boundary_blocks)  v.clear();
                
                const auto create_geom = [&](const std::size_t& lim, const auto& loc_glob, auto& geom)
                {
                    for (std::size_t lb = 0; lb < lim; ++lb)
                    {
                        const auto bidx           = utils::tag[loc_glob](lb);
                        const std::size_t lb_glob = grid_partition.to_global(bidx).value;
                        const auto box            = block_arrangement.get_block_box(lb_glob);
                        const auto ibdy           = block_arrangement.is_domain_boundary(lb_glob);
                        geom.bounding_boxes.push_back(box);
                        geom.domain_boundary.push_back(ibdy);
                        
                        for (int i = 0; i < 6; ++i)
                        {
                            if (ibdy.bnds[i])
                            {
                                geom.boundary_blocks[i].push_back(lb);
                            }
                        }
                        
                    }
                };
                
                create_geom(grid_partition.get_num_global_blocks(), partition::global, global_geometry);
                create_geom(grid_partition.get_num_local_blocks(),  partition::local,  local_geometry);
                
                global_geometry.bounding_boxes.transfer();
                global_geometry.domain_boundary.transfer();
                for (int i = 0; i < 6; ++i) global_geometry.boundary_blocks[i].transfer();
                
                local_geometry.bounding_boxes.transfer();
                local_geometry.domain_boundary.transfer();
                for (int i = 0; i < 6; ++i) local_geometry.boundary_blocks[i].transfer();
            }
            
            template <typename loc_glob_t, typename device_t>
            const geomety_image_type image(const loc_glob_t& i, const device_t& dev) const
            {
                const auto& geom = this->geometry(i);
                geomety_image_type output(geom.coords, geom.num_cell);
                output.bounding_boxes  = utils::make_vec_image(geom.bounding_boxes.data(dev));
                output.domain_boundary = utils::make_vec_image(geom.domain_boundary.data(dev));
                for (int i = 0; i < 6; ++i)
                {
                    output.boundary_blocks[i] = utils::make_vec_image(geom.boundary_blocks[i].data(dev));
                }
                return output;
            }
            
            template <typename device_t>
            const geomety_image_type image(const device_t& dev) const
            {
                return this->image(partition::local, dev);
            }
            
            constexpr static int dim() {return grid_dim;}
            
            std::size_t  get_grid_size()                                const { return get_num_cells(0)*get_num_cells(1)*get_num_cells(2)*get_num_global_blocks(); }
            std::size_t  get_num_local_blocks()                         const { return grid_partition.get_num_local_blocks(); }
            std::size_t  get_num_global_blocks()                        const { return grid_partition.get_num_global_blocks(); }
            std::size_t  get_num_blocks(const partition::global_tag_t&) const { return grid_partition.get_num_global_blocks(); }
            std::size_t  get_num_blocks(const partition::local_tag_t&)  const { return grid_partition.get_num_local_blocks(); }
            int get_num_cells(const std::size_t& i)                     const { return global_geometry.num_cell[i]; }
            bound_box_t<dtype,  3> get_bounds()                         const { return block_arrangement.get_bounds(); }
            const par_group_t& group()                                  const { return grid_group; }
            const coord_t& coord_sys()                                  const { return global_geometry.get_coords(); }
            const auto& get_blocks()                                    const { return block_arrangement; }
            auto& get_blocks()                                                { return block_arrangement; }
            const partition::block_partition_t& get_partition()         const { return grid_partition; }
            constexpr static bool is_3d()                                     { return dim()==3; }
            
            template <typename idx_t> _finline_ coord_point_type get_coords(const idx_t& i) const
            {
                return global_geometry.get_coords(i);
            }
                
            template <typename idx_t> _finline_ coord_point_type get_comp_coords(const idx_t& i) const
            {
                return global_geometry.get_comp_coords(i);
            }
            
            const auto is_domain_boundary(const partition::partition_tagged auto& lb) const
            {
                return global_geometry.is_domain_boundary(grid_partition.to_global(lb).value);
            }
            
            ctrs::array<dtype,  3> get_dx(const std::size_t& lb) const
            { 
                return this->get_dx(utils::tag[partition::local](lb));
            }
            
            dtype get_dx(const int i, const std::size_t& lb)  const
            {
                return this->get_dx(i, utils::tag[partition::local](lb));
            }
            
            template <partition::partition_tagged idx_t>
            ctrs::array<dtype,  3> get_dx(const idx_t& lb) const
            {
                return this->geometry(typename idx_t::tag_type()).get_dx(lb.value);
            }
            
            template <partition::partition_tagged idx_t>
            dtype get_dx(const int i, const idx_t& lb)  const
            {
                return this->geometry(typename idx_t::tag_type()).get_dx(i, lb.value);
            }
            
            template <partition::partition_tagged idx_t>
            bound_box_t<dtype, 3> get_bounding_box(const idx_t& lb) const
            {
                return this->geometry(typename idx_t::tag_type()).get_bounding_box(grid_partition.to_global(lb).value);
            }
            
            bound_box_t<dtype, 3> get_bounding_box(const std::size_t& lb) const
            {
                return this->get_bounding_box(utils::tag[partition::global](lb));
            }
            
            auto compute_dx_min() const
            {
                ctrs::array<dtype,  3> output(1e50, 1e50, 1e50);
                for (std::size_t lb = 0; lb < this->get_num_global_blocks(); ++lb)
                {
                    auto dx = this->get_dx(lb);
                    for (int i = 0; i < output.size(); ++i)
                    {
                        output[i] = utils::min(output[i], dx[i]);
                    }
                }
                return output;
            }
            
            const auto& get_coord_sys() const {return global_geometry.get_coords();}
            
            const auto& geometry(const partition::global_tag_t&) const { return global_geometry; }
            const auto& geometry(const partition::local_tag_t&)  const { return local_geometry;  }
            
            template <typename container_t, typename func_t>
            requires (partition::partition_tagged<typename container_t::value_type>)
            void select_blocks(container_t& output, const func_t& func) const
            {
                output.clear();
                constexpr auto log = []
                {
                    if constexpr (std::same_as<typename container_t::value_type::tag_type, partition::local_tag_t>) return partition::local;
                    else return partition::global;
                }();
                
                for (std::size_t lb = 0; lb < get_num_blocks(log); ++lb)
                {
                    auto tag_lb = utils::tag[log](lb);
                    if (func(tag_lb)) output.push_back(tag_lb);
                }
                
                if constexpr ( requires { output.transfer(); } ) output.transfer();
            }
            
            template <typename container_t, typename func_t, typename loc_or_glob_t>
            void select_blocks(container_t& output, const func_t& func, const loc_or_glob_t& log) const
            {
                output.clear();
                for (std::size_t lb = 0; lb < get_num_blocks(log); ++lb)
                {
                    auto tag_lb = utils::tag[log](lb);
                    if (func(tag_lb)) output.push_back(tag_lb);
                }
                if constexpr ( requires { output.transfer(); } ) output.transfer();
            }
            
            template <typename func_t, typename loc_or_glob_t>
            auto select_blocks(const func_t& func, const loc_or_glob_t& log) const
            {                
                using idx_t = decltype(utils::tag[log](std::size_t(0)));
                std::vector<idx_t> output;
                this->select_blocks(output, func);
                return output;
            }
            
            static constexpr bool periodic_refinement_default = false;
            
            template <typename list_t, typename refs_t, typename constraint_t>
            void refine_blocks(const list_t& list, ctrs::array<bool, dim()>& periodic, const refs_t& refs, const constraint_t& constraint = amr::constraints::factor2)
            {
                std::size_t vsize = 1;
                constexpr bool list_is_vec = utils::is_std_vector<list_t>;
                constexpr bool refs_is_vec = utils::is_std_vector<refs_t>;
                if constexpr (list_is_vec) vsize = utils::max(vsize, list.size());
                if constexpr (refs_is_vec) vsize = utils::max(vsize, refs.size());
                auto arg_list = [&]()
                {
                    if constexpr (list_is_vec) return list;
                    else return std::vector<list_t>(vsize, list);
                }();
                
                const auto ref_list = [&]()
                {
                    if constexpr (refs_is_vec) return refs;
                    else return std::vector<refs_t>(vsize, refs);
                }();
                
                using idx_t = decltype(arg_list)::value_type;
                if constexpr (partition::partition_tagged<idx_t>)
                {
                    static_assert(std::same_as<typename idx_t::tag_type, partition::global_tag_t>, "Not allowed to refine on a locally-partitioned basis");
                    std::vector<std::size_t> rdata;
                    rdata.reserve(list.size());
                    //Need to convert the tagged integer type to a basic integer type
                    for (const auto& plb: list) rdata.push_back(plb);
                    this->get_blocks().refine(rdata, periodic, ref_list, constraint);
                }
                else
                {
                    this->get_blocks().refine(list, periodic, ref_list, constraint);
                }
                compute_geometry();
                
                for (auto& p: dependents) p->on_blocks_update();
            }
            
            template <typename list_t, typename constraint_t = decltype(amr::constraints::factor2)>
            void refine_blocks(const list_t& list, ctrs::array<bool, dim()>& periodic, const constraint_t& constraint = amr::constraints::factor2)
            {
                const typename block_arrangement_t::refine_type rtype = true;
                this->refine_blocks(list, periodic, rtype, constraint);
            }
            
            template <typename list_t>
            void refine_blocks(const list_t& list)
            {
                const typename block_arrangement_t::refine_type rtype = true;
                ctrs::array<bool, dim()> periodic = periodic_refinement_default;
                this->refine_blocks(list, periodic, rtype, amr::constraints::factor2);
            }
    };
}