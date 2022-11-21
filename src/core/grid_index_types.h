#pragma once
#include <type_traits>

#include "core/config.h"
#include "core/typedef.h"
#include "core/ctrs.h"
namespace spade::grid
{
    enum array_centering
    {
        cell_centered,
        node_centered,
        face_centered,
        edge_centered
    };
    

    struct cell_idx_t : public ctrs::arithmetic_array_t<int, 4, cell_idx_t>
    {
        constexpr static array_centering centering_type() {return cell_centered;}
        using base_t = ctrs::arithmetic_array_t<int, 4, cell_idx_t>;
        using base_t::base_t;
        
        static constexpr int i_idx =   0;
        static constexpr int j_idx =   1;
        static constexpr int k_idx =   2;
        static constexpr int lb_idx =  3;
        
        cell_idx_t(){}
        int& i () {return (*this)[i_idx];}
        int& i (const int& idx) {return (*this)[i_idx + idx];}
        int& j () {return (*this)[j_idx];}
        int& k () {return (*this)[k_idx];}
        int& lb() {return (*this)[lb_idx];}
        const int& i () const {return (*this)[i_idx];}
        const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        const int& j () const {return (*this)[j_idx];}
        const int& k () const {return (*this)[k_idx];}
        const int& lb() const {return (*this)[lb_idx];}
    };
    
    struct face_idx_t : public ctrs::arithmetic_array_t<int, 5, face_idx_t>
    {
        constexpr static array_centering centering_type() {return face_centered;}
        using base_t = ctrs::arithmetic_array_t<int, 5, face_idx_t>;
        using base_t::base_t;
        
        static constexpr int dir_idx = 0;
        static constexpr int i_idx =   1;
        static constexpr int j_idx =   2;
        static constexpr int k_idx =   3;
        static constexpr int lb_idx =  4;
        
        face_idx_t(){}
        int& dir() {return (*this)[dir_idx];}
        int& i  () {return (*this)[i_idx];}
        int& i (const int& idx) {return (*this)[i_idx + idx];}
        int& j  () {return (*this)[j_idx];}
        int& k  () {return (*this)[k_idx];}
        int& lb () {return (*this)[lb_idx];}
        const int& dir() const {return (*this)[dir_idx];}
        const int& i  () const {return (*this)[i_idx];}
        const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        const int& j  () const {return (*this)[j_idx];}
        const int& k  () const {return (*this)[k_idx];}
        const int& lb () const {return (*this)[lb_idx];}
    };
    
    struct node_idx_t : public ctrs::arithmetic_array_t<int, 4, node_idx_t>
    {
        constexpr static array_centering centering_type() {return node_centered;}
        using base_t = ctrs::arithmetic_array_t<int, 4, node_idx_t>;
        using base_t::base_t;
        
        static constexpr int i_idx =   0;
        static constexpr int j_idx =   1;
        static constexpr int k_idx =   2;
        static constexpr int lb_idx =  3;
        
        node_idx_t(){}
        int& i () {return (*this)[i_idx];}
        int& i (const int& idx) {return (*this)[i_idx + idx];}
        int& j () {return (*this)[j_idx];}
        int& k () {return (*this)[k_idx];}
        int& lb() {return (*this)[lb_idx];}
        const int& i () const {return (*this)[i_idx];}
        const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        const int& j () const {return (*this)[j_idx];}
        const int& k () const {return (*this)[k_idx];}
        const int& lb() const {return (*this)[lb_idx];}
    };
    
    //i_cell is the cell index
    //idir is the principal cartesian component of the face normal vector [unsigned]
    //pm = 0 indicates the negative-side face, pm = 1 indicates the positive-side face
    static face_idx_t cell_to_face(const cell_idx_t& i_cell, const int& idir, const int& pm)
    {
        face_idx_t output;
        output.i()   = i_cell.i();
        output.j()   = i_cell.j();
        output.k()   = i_cell.k();
        output.lb()  = i_cell.lb();
        output.dir() = idir;
        
        //May cause an issue if i,j,k are not stored in consecution
        output.i(idir) += pm;
        return output;
    }
    
    //i_cell is the cell index
    //offset has 3 entries that are either 0 or 1, indicating which corner of the cell to be chosen
    static node_idx_t cell_to_node(const cell_idx_t& i_cell, const ctrs::array<int,3>& offset)
    {
        print("NOT IMPLEMENTED:", __FILE__, __LINE__);
        abort();
        // node_idx_t output((int)i_cell[0], (int)i_cell[1], (int)i_cell[2], (int)i_cell[3]);
        node_idx_t output;
        // output[0] += offset[0];
        // output[1] += offset[1];
        // output[2] += offset[2];
        return output;
    }
    
    //i_face is the face index
    //pm = 0 indicates the negative-side ("left") face, pm = 1 indicates the positive-side ("right") face
    static cell_idx_t face_to_cell(const face_idx_t& i_face, const int& pm)
    {
        cell_idx_t output;
        output.i() = i_face.i();
        output.j() = i_face.j();
        output.k() = i_face.k();
        output.lb() = i_face.lb();
        output.i(i_face.dir()) += (pm-1);

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
            (rtype)i_cell.i()+0.5,
            (rtype)i_cell.j()+0.5,
            (rtype)i_cell.k()+0.5);
        return output;
    }
    
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const node_idx_t& i_node)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_node.i(),
            (rtype)i_node.j(),
            (rtype)i_node.k());
        return output;
    }
    
    template <typename rtype=real_t> ctrs::array<rtype, 3> get_index_coord(const face_idx_t& i_face)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_face.i()+0.5,
            (rtype)i_face.j()+0.5,
            (rtype)i_face.k()+0.5);
        output[i_face.dir_idx] -= 0.5;
        return output;
    }
}