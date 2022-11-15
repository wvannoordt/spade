#pragma once
#include <type_traits>
#include <concepts>

#include "core/typedef.h"
#include "core/ctrs.h"
namespace spade::grid
{
    enum array_center_e
    {
        cell_centered,
        node_centered,
        face_centered,
        edge_centered
    };
    
    using index_integral_base_t = int;
    
    template <typename idx_t = index_integral_base_t>
    struct cell_idx_t : public ctrs::arithmetic_array_t<idx_t, 4, cell_idx_t<idx_t>>
    {
        constexpr static array_center_e centering_type() {return cell_centered;}
        using base_t = ctrs::arithmetic_array_t<idx_t, 4, cell_idx_t<idx_t>>;
        using base_t::base_t;
        cell_idx_t(){}
        idx_t& i () {return (*this)[0];}
        idx_t& j () {return (*this)[1];}
        idx_t& k () {return (*this)[2];}
        idx_t& lb() {return (*this)[3];}
        const idx_t& i () const {return (*this)[0];}
        const idx_t& j () const {return (*this)[1];}
        const idx_t& k () const {return (*this)[2];}
        const idx_t& lb() const {return (*this)[3];}
    };
    
    template <typename idx_t = index_integral_base_t>
    struct face_idx_t : public ctrs::arithmetic_array_t<idx_t, 5, face_idx_t<idx_t>>
    {
        constexpr static array_center_e centering_type() {return face_centered;}
        using base_t = ctrs::arithmetic_array_t<idx_t, 5, face_idx_t<idx_t>>;
        using base_t::base_t;
        face_idx_t(){}
        idx_t& dir() {return (*this)[0];}
        idx_t& i  () {return (*this)[1];}
        idx_t& j  () {return (*this)[2];}
        idx_t& k  () {return (*this)[3];}
        idx_t& lb () {return (*this)[4];}
        const idx_t& dir() const {return (*this)[0];}
        const idx_t& i  () const {return (*this)[1];}
        const idx_t& j  () const {return (*this)[2];}
        const idx_t& k  () const {return (*this)[3];}
        const idx_t& lb () const {return (*this)[4];}
    };
    
    template <typename idx_t = index_integral_base_t>
    struct node_idx_t : public ctrs::arithmetic_array_t<idx_t, 4, node_idx_t<idx_t>>
    {
        constexpr static array_center_e centering_type() {return node_centered;}
        using base_t = ctrs::arithmetic_array_t<idx_t, 4, node_idx_t<idx_t>>;
        using base_t::base_t;
        node_idx_t(){}
        idx_t& i  () {return (*this)[0];}
        idx_t& j  () {return (*this)[1];}
        idx_t& k  () {return (*this)[2];}
        idx_t& lb () {return (*this)[3];}
        const idx_t& i  () const {return (*this)[0];}
        const idx_t& j  () const {return (*this)[1];}
        const idx_t& k  () const {return (*this)[2];}
        const idx_t& lb () const {return (*this)[3];}
    };
    
    template <typename T> concept multiblock_grid_idx_t = requires(T t)
    {
        t;
    };
    
    enum subindex_e
    {
        i_subindex,
        j_subindex,
        k_subindex,
        lb_subindex,
        dir_subindex
    };
    
    template <const array_center_e centering, const subindex_e subidx> struct get_subidx
    {
        //specializations only
    };
    
    //NOTE: it is currently required that i,j,k are stored next to each other!
    template<> struct get_subidx<cell_centered, i_subindex>   {static constexpr int value =  0;};
    template<> struct get_subidx<cell_centered, j_subindex>   {static constexpr int value =  1;};
    template<> struct get_subidx<cell_centered, k_subindex>   {static constexpr int value =  2;};
    template<> struct get_subidx<cell_centered, lb_subindex>  {static constexpr int value =  3;};
    template<> struct get_subidx<cell_centered, dir_subindex> {}; //invalid
    
    template<> struct get_subidx<face_centered, i_subindex>   {static constexpr int value =  1;};
    template<> struct get_subidx<face_centered, j_subindex>   {static constexpr int value =  2;};
    template<> struct get_subidx<face_centered, k_subindex>   {static constexpr int value =  3;};
    template<> struct get_subidx<face_centered, lb_subindex>  {static constexpr int value =  4;};
    template<> struct get_subidx<face_centered, dir_subindex> {static constexpr int value =  0;};
    
    template<> struct get_subidx<node_centered, i_subindex>   {static constexpr int value =  0;};
    template<> struct get_subidx<node_centered, j_subindex>   {static constexpr int value =  1;};
    template<> struct get_subidx<node_centered, k_subindex>   {static constexpr int value =  2;};
    template<> struct get_subidx<node_centered, lb_subindex>  {static constexpr int value =  3;};
    template<> struct get_subidx<node_centered, dir_subindex> {}; //invalid
    
    static constexpr int get_face_dir_idx()
    {
        return get_subidx<face_centered, dir_subindex>::value;
    }
    
    static int get_face_dir(const face_idx_t& i_face)
    {
        return i_face[get_subidx<face_centered, dir_subindex>::value];
    }
    
    
    static int get_block_number(const cell_idx_t& i_cell)
    {
        return i_cell[get_subidx<cell_centered, lb_subindex>::value];
    }
    
    static int get_block_number(const face_idx_t& i_face)
    {
        return i_face[get_subidx<face_centered, lb_subindex>::value];
    }
    
    static int get_block_number(const node_idx_t& i_node)
    {
        return i_node[get_subidx<node_centered, lb_subindex>::value];
    }
    
    //i_cell is the cell index
    //idir is the principal cartesian component of the face normal vector [unsigned]
    //pm = 0 indicates the negative-side face, pm = 1 indicates the positive-side face
    static face_idx_t cell_to_face(const cell_idx_t& i_cell, const int& idir, const int& pm)
    {
        face_idx_t output;
        output[get_subidx<face_centered,   i_subindex>::value] = (int)i_cell[get_subidx<cell_centered,  i_subindex>::value];
        output[get_subidx<face_centered,   j_subindex>::value] = (int)i_cell[get_subidx<cell_centered,  j_subindex>::value];
        output[get_subidx<face_centered,   k_subindex>::value] = (int)i_cell[get_subidx<cell_centered,  k_subindex>::value];
        output[get_subidx<face_centered,  lb_subindex>::value] = (int)i_cell[get_subidx<cell_centered, lb_subindex>::value];
        output[get_subidx<face_centered, dir_subindex>::value] = idir;
        
        //May cause an issue if i,j,k are not stored in consecution
        output[get_subidx<face_centered,  i_subindex>::value + idir] += pm;
        return output;
    }
    
    //i_cell is the cell index
    //offset has 3 entries that are either 0 or 1, indicating which corner of the cell to be chosen
    static node_idx_t cell_to_node(const cell_idx_t& i_cell, const ctrs::array<int,3>& offset)
    {
        print("NOT IMPLEMENTED:", __FILE__, __LINE__);
        abort();
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
        cell_idx_t output;
        output[get_subidx<cell_centered,  i_subindex>::value] = (int)i_face[get_subidx<face_centered,  i_subindex>::value];
        output[get_subidx<cell_centered,  j_subindex>::value] = (int)i_face[get_subidx<face_centered,  j_subindex>::value];
        output[get_subidx<cell_centered,  k_subindex>::value] = (int)i_face[get_subidx<face_centered,  k_subindex>::value];
        output[get_subidx<cell_centered, lb_subindex>::value] = (int)i_face[get_subidx<face_centered, lb_subindex>::value];
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
        ctrs::array<rtype, 3> output(
            (rtype)i_cell[get_subidx<cell_centered, i_subindex>::value]+0.5,
            (rtype)i_cell[get_subidx<cell_centered, j_subindex>::value]+0.5,
            (rtype)i_cell[get_subidx<cell_centered, k_subindex>::value]+0.5);
        return output;
    }
    
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const node_idx_t& i_node)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_node[get_subidx<node_centered, i_subindex>::value],
            (rtype)i_node[get_subidx<node_centered, j_subindex>::value],
            (rtype)i_node[get_subidx<node_centered, k_subindex>::value]);
        return output;
    }
    
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const face_idx_t& i_face)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_face[get_subidx<face_centered, i_subindex>::value]+0.5,
            (rtype)i_face[get_subidx<face_centered, j_subindex>::value]+0.5,
            (rtype)i_face[get_subidx<face_centered, k_subindex>::value]+0.5);
        output[get_face_dir(i_face)] -= 0.5;
        return output;
    }
    
    template <const array_center_e centering>
    requires (centering==cell_centered)
    static cell_idx_t make_idx(const int& i, const int& j, const int& k, const int& lb)
    {
        cell_idx_t output;
        output[get_subidx<cell_centered,  i_subindex>::value] = i;
        output[get_subidx<cell_centered,  j_subindex>::value] = j;
        output[get_subidx<cell_centered,  k_subindex>::value] = k;
        output[get_subidx<cell_centered, lb_subindex>::value] = lb;
        return output;
    }
    
    template <const array_center_e centering>
    requires (centering==face_centered)
    static face_idx_t make_idx(const int& idir, const int& i, const int& j, const int& k, const int& lb)
    {
        face_idx_t output;
        output[get_subidx<face_centered,dir_subindex>::value] = idir;
        output[get_subidx<face_centered,  i_subindex>::value] = i;
        output[get_subidx<face_centered,  j_subindex>::value] = j;
        output[get_subidx<face_centered,  k_subindex>::value] = k;
        output[get_subidx<face_centered, lb_subindex>::value] = lb;
        return output;
    }
    
    template <const array_center_e centering>
    requires (centering==node_centered)
    static node_idx_t make_idx(const int& i, const int& j, const int& k, const int& lb)
    {
        node_idx_t output;
        output[get_subidx<node_centered,  i_subindex>::value] = i;
        output[get_subidx<node_centered,  j_subindex>::value] = j;
        output[get_subidx<node_centered,  k_subindex>::value] = k;
        output[get_subidx<node_centered, lb_subindex>::value] = lb;
        return output;
    }
}