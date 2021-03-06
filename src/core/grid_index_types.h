#pragma once
#include <type_traits>
#include <concepts>

#include "core/typedef.h"
#include "core/ctrs.h"
namespace spade::grid
{
    
    template <typename derived_t, typename arith_t>
    struct arithmetic_wrapper_t
    {
        arith_t val;
        arithmetic_wrapper_t(void) {}
        arithmetic_wrapper_t(const arithmetic_wrapper_t& rhs) {val = rhs.val;}
        arithmetic_wrapper_t(const arith_t& rhs) {val = rhs;}
        
        derived_t& operator = (const arith_t& rhs)
        {
            val = rhs;
            return static_cast<derived_t&>(*this);
        }
        
        derived_t& operator /= (const arith_t& rhs)
        {
            val /= rhs;
            return static_cast<derived_t&>(*this);
        }
        derived_t& operator *= (const arith_t& rhs)
        {
            val *= rhs;
            return static_cast<derived_t&>(*this);
        }
        
        derived_t& operator += (const arith_t& rhs)
        {
            val += rhs;
            return static_cast<derived_t&>(*this);
        }
        derived_t& operator -= (const arith_t& rhs)
        {
            val -= rhs;
            return static_cast<derived_t&>(*this);
        }
        operator arith_t() const { return val; }
    };
    
    template <typename arith_t> struct cell_t : public arithmetic_wrapper_t<cell_t<arith_t>, arith_t>
    {
        typedef arithmetic_wrapper_t<cell_t<arith_t>, arith_t> base_t;
        typedef arith_t value_type;
        using base_t::base_t;
    };
    template <typename arith_t> struct node_t : public arithmetic_wrapper_t<node_t<arith_t>, arith_t>
    {
        typedef arithmetic_wrapper_t<node_t<arith_t>, arith_t> base_t;
        typedef arith_t value_type;
        using base_t::base_t;
    };
    template <typename arith_t> struct face_t : public arithmetic_wrapper_t<face_t<arith_t>, arith_t>
    {
        typedef arithmetic_wrapper_t<face_t<arith_t>, arith_t> base_t;
        typedef arith_t value_type;
        using base_t::base_t;
    };
    
    typedef ctrs::array<cell_t<int>, 4> cell_idx_t;
    typedef ctrs::array<face_t<int>, 5> face_idx_t;
    typedef ctrs::array<node_t<int>, 4> node_idx_t;
    
    template <typename T> concept multiblock_grid_idx_t = requires(T t)
    {
        t;
    };
    
    static auto get_face_dir_idx(void)
    {
        return 4;
    }
    
    static auto get_face_dir(const face_idx_t& i_face)
    {
        return i_face[get_face_dir_idx()];
    }
    
    //i_cell is the cell index
    //idir is the principal cartesian component of the face normal vector [unsigned]
    //pm = 0 indicates the negative-side face, pm = 1 indicates the positive-side face
    static face_idx_t cell_to_face(const cell_idx_t& i_cell, const int& idir, const int& pm)
    {
        //Note that the direction is presently the last index
        face_idx_t output((int)i_cell[0], (int)i_cell[1], (int)i_cell[2], (int)i_cell[3], idir);
        output[idir] += pm;
        return output;
    }
    
    //i_cell is the cell index
    //offset has 3 entries that are either 0 or 1, indicating which corner of the cell to be chosen
    static node_idx_t cell_to_node(const cell_idx_t& i_cell, const ctrs::array<int,3>& offset)
    {
        node_idx_t output((int)i_cell[0], (int)i_cell[1], (int)i_cell[2], (int)i_cell[3]);
        output[0] += offset[0];
        output[1] += offset[1];
        output[2] += offset[2];
        return output;
    }
    
    //i_face is the face index
    //pm = 0 indicates the negative-side ("left") face, pm = 1 indicates the positive-side ("right") face
    static cell_idx_t face_to_cell(const face_idx_t& i_face, const int& pm)
    {
        cell_idx_t output((int)i_face[0], (int)i_face[1], (int)i_face[2], (int)i_face[3]);
        output[get_face_dir(i_face)] += (pm-1);
        return output;
    }
    
    //i_face is the face index
    //pm = 0 indicates the negative-side ("left") face, pm = 1 indicates the positive-side ("right") face
    static node_idx_t face_to_node(const face_idx_t& i_face, const int& pm)
    {
        print("NOT IMPLEMENTED:", __FILE__, __LINE__);
        abort();
        node_idx_t output((int)i_face[0], (int)i_face[1], (int)i_face[2], (int)i_face[3]);
        return output;
    }
    
    //Return the index-space coordinates of the index, used to compute the computational coordinates
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const cell_idx_t& i_cell)
    {
        ctrs::array<rtype, 3> output((rtype)i_cell[0]+0.5, (rtype)i_cell[1]+0.5, (rtype)i_cell[2]+0.5);
        return output;
    }
    
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const node_idx_t& i_node)
    {
        ctrs::array<rtype, 3> output((rtype)i_node[0], (rtype)i_node[1], (rtype)i_node[2]);
        return output;
    }
    
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const face_idx_t& i_face)
    {
        ctrs::array<rtype, 3> output((rtype)i_face[0]+0.5, (rtype)i_face[1]+0.5, (rtype)i_face[2]+0.5);
        output[get_face_dir(i_face)] -= 0.5;
        return output;
    }
}