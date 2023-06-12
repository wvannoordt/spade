#pragma once

#include "grid/grid_index_types.h"
#include "core/static_math.h"

namespace spade::omni
{
    namespace detail
    {
        //Yuck, there is a lot of garbage here
        template <const int ii_seek, const int ii_spec, const grid::array_centering ctr, const grid::array_centering... ctrs>
        struct get_val_helper_t
        {
            constexpr static grid::array_centering value = get_val_helper_t<ii_seek+1, ii_spec, ctrs...>::value;
        };
        
        template <const int ii, const grid::array_centering ctr, const grid::array_centering... ctrs>
        struct get_val_helper_t<ii, ii, ctr, ctrs...>
        {
            constexpr static grid::array_centering value = ctr;
        };
        
        template <const int ii, const grid::array_centering... ctrs>
        struct get_val_t
        {
            constexpr static grid::array_centering value = get_val_helper_t<0, ii, ctrs...>::value;
        };
        
        template <const grid::array_centering... ctr>
        struct ctable_t
        {
            template <const int ii> constexpr static grid::array_centering elem = get_val_t<ii, ctr...>::value;
            constexpr static int size() {return sizeof...(ctr);}
        };
        
        using ctable_inst_t = ctable_t
        <
            grid::cell_centered,
            grid::face_centered,
            grid::edge_centered,
            grid::node_centered,
            grid::edge_centered,
            grid::face_centered
        >;
        
        template <const grid::array_centering ctr, const int mod_di, const int mod_dj, const int mod_dk>
        struct compute_center_t
        {
            constexpr static int i_base = ctr==grid::node_centered?3:0;
            constexpr static grid::array_centering centering()
            {
                //just trust the algorithm
                return ctable_inst_t::elem<static_math::mod<i_base + mod_di + mod_dj + mod_dk, ctable_inst_t::size()>::value>;
            }
        };
        
        //Note that we are guaranteed to have a 1-d stencil here.
        template <const int mod_di, const int mod_dj, const int mod_dk> struct compute_center_t<grid::face_centered, mod_di, mod_dj, mod_dk>
        {
            constexpr static grid::array_centering centering() {return (mod_di == 1)?grid::cell_centered:grid::face_centered;}
        };
    }
    
    template <const int di, const int dj, const int dk> struct offset_t
    {
        constexpr static bool is_1d = ((dj==0) && (dk==0));
        
        //if we are at a node of type 'center', what kind of node is at the given offset?
        template <const grid::array_centering center>
        constexpr static grid::array_centering relative_node_centering
            = detail::compute_center_t<
                center,
                static_math::mod<di, 2>::value,
                static_math::mod<dj, 2>::value,
                static_math::mod<dk, 2>::value
                >::centering();

        template<grid::grid_index index_t>
        requires((index_t::centering_type() == grid::cell_centered) || (index_t::centering_type() == grid::node_centered))
        constexpr static auto compute_index(const index_t& i)
        {
            const int di_div = static_math::moddiv<di,2>::value;
            const int dj_div = static_math::moddiv<dj,2>::value;
            const int dk_div = static_math::moddiv<dk,2>::value;
            
            const int di_mod = static_math::mod<di,2>::value;
            const int dj_mod = static_math::mod<dj,2>::value;
            const int dk_mod = static_math::mod<dk,2>::value;
            
            using output_type = typename grid::get_index_type<relative_node_centering<index_t::centering_type()>>::array_type;
            static_assert(output_type::centering_type() != grid::edge_centered, "edge centering not implemented");
            static_assert(output_type::centering_type() != grid::node_centered, "node centering not implemented");
            // START HERE
            index_t output;
            output.lb() = i.lb();
            output.i() = i.i() + di_div;
            output.j() = i.j() + dj_div;
            output.k() = i.k() + dk_div;
            if constexpr (output_type::centering_type() == grid::face_centered)
            {
                return cell_to_face(output, dj_mod+2*dk_mod, 1);
            }
            else if constexpr (output_type::centering_type() == grid::node_centered)
            {
                //GARBAGE! needs implementing
                return cell_to_node(output,ctrs::array<int,3>());
            }
            else
            {
                return output;
            }
        }
        
        template<grid::grid_index index_t>
        requires(((index_t::centering_type() == grid::face_centered) || (index_t::centering_type() == grid::edge_centered)) && is_1d)
        constexpr static auto compute_index(const index_t& i)
        {
            const int di_div = static_math::moddiv<di,2>::value;
            const int di_mod = static_math::mod<di,2>::value;
            
            static_assert(index_t::centering_type() != grid::edge_centered, "edge centering not implemented");
            
            using output_type = typename grid::get_index_type<relative_node_centering<index_t::centering_type()>>::array_type;
            static_assert(output_type::centering_type() != grid::edge_centered, "edge centering not implemented");
            index_t output;
            output.lb()    = i.lb();
            output.i()     = i.i();
            output.j()     = i.j();
            output.k()     = i.k();
            output.dir()   = i.dir();
            output.i(i.dir()) += di_div;
            if constexpr (output_type::centering_type() == grid::cell_centered)
            {
                return face_to_cell(output, 1);
            }
            else
            {
                return output;
            }
        }
    };
}