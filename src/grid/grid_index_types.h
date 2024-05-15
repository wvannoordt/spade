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
        edge_centered,
        agno_centered  //'agnostic' centering
    };
    
    template <const array_centering ctr> struct get_complement_centering {};
    template <> struct get_complement_centering<edge_centered> {constexpr static array_centering value = node_centered;};
    template <> struct get_complement_centering<face_centered> {constexpr static array_centering value = cell_centered;};
    template <> struct get_complement_centering<node_centered> {constexpr static array_centering value = edge_centered;};
    template <> struct get_complement_centering<cell_centered> {constexpr static array_centering value = face_centered;};
    
    template <const array_centering ctr> struct get_dual_centering {};
    template <> struct get_dual_centering<edge_centered> {constexpr static array_centering value = face_centered;};
    template <> struct get_dual_centering<face_centered> {constexpr static array_centering value = edge_centered;};
    template <> struct get_dual_centering<node_centered> {constexpr static array_centering value = cell_centered;};
    template <> struct get_dual_centering<cell_centered> {constexpr static array_centering value = node_centered;};
    

    template <typename T> concept grid_index = requires(T t)
    {
        { T::centering_type() } -> std::same_as<array_centering>;
    };
    
    template <typename T> concept is_directed_index = grid_index<T> && ( T::centering_type() == face_centered || T::centering_type() == edge_centered);
    
    struct cell_idx_t : public ctrs::arithmetic_array_t<int, 4, cell_idx_t>
    {
        constexpr static array_centering centering_type() {return cell_centered;}
        using base_t = ctrs::arithmetic_array_t<int, 4, cell_idx_t>;
        using base_t::base_t;
        
        static constexpr int i_idx =   0;
        static constexpr int j_idx =   1;
        static constexpr int k_idx =   2;
        static constexpr int lb_idx =  3;
        
        _sp_hybrid cell_idx_t(){}
        _sp_hybrid int& i () {return (*this)[i_idx];}
        _sp_hybrid int& i (const int& idx) {return (*this)[i_idx + idx];}
        _sp_hybrid int& j () {return (*this)[j_idx];}
        _sp_hybrid int& k () {return (*this)[k_idx];}
        _sp_hybrid int& lb() {return (*this)[lb_idx];}
        _sp_hybrid const int& i () const {return (*this)[i_idx];}
        _sp_hybrid const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        _sp_hybrid const int& j () const {return (*this)[j_idx];}
        _sp_hybrid const int& k () const {return (*this)[k_idx];}
        _sp_hybrid const int& lb() const {return (*this)[lb_idx];}

        bool operator == (const cell_idx_t& rhs) const
        {
            return (i() == rhs.i()) && (j() == rhs.j()) && (k() == rhs.k()) && (lb() == rhs.lb());
        }
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
        
        _sp_hybrid face_idx_t(){}
        _sp_hybrid int& dir() {return (*this)[dir_idx];}
        _sp_hybrid int& i  () {return (*this)[i_idx];}
        _sp_hybrid int& i (const int& idx) {return (*this)[i_idx + idx];}
        _sp_hybrid int& j  () {return (*this)[j_idx];}
        _sp_hybrid int& k  () {return (*this)[k_idx];}
        _sp_hybrid int& lb () {return (*this)[lb_idx];}
        _sp_hybrid const int& dir() const {return (*this)[dir_idx];}
        _sp_hybrid const int& i  () const {return (*this)[i_idx];}
        _sp_hybrid const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        _sp_hybrid const int& j  () const {return (*this)[j_idx];}
        _sp_hybrid const int& k  () const {return (*this)[k_idx];}
        _sp_hybrid const int& lb () const {return (*this)[lb_idx];}

        bool operator == (const face_idx_t& rhs) const
        {
            return (i() == rhs.i()) && (j() == rhs.j()) && (k() == rhs.k()) && (lb() == rhs.lb()) && (dir() == rhs.dir());
        }
    };
    
    struct edge_idx_t : public ctrs::arithmetic_array_t<int, 5, edge_idx_t>
    {
        //CURRENTLY UNIMPLEMENTED!!!
        constexpr static array_centering centering_type() {return edge_centered;}
        using base_t = ctrs::arithmetic_array_t<int, 5, edge_idx_t>;
        using base_t::base_t;
        
        static constexpr int dir_idx = 0;
        static constexpr int i_idx =   1;
        static constexpr int j_idx =   2;
        static constexpr int k_idx =   3;
        static constexpr int lb_idx =  4;
        
        _sp_hybrid edge_idx_t(){}
        _sp_hybrid int& dir() {return (*this)[dir_idx];}
        _sp_hybrid int& i  () {return (*this)[i_idx];}
        _sp_hybrid int& i (const int& idx) {return (*this)[i_idx + idx];}
        _sp_hybrid int& j  () {return (*this)[j_idx];}
        _sp_hybrid int& k  () {return (*this)[k_idx];}
        _sp_hybrid int& lb () {return (*this)[lb_idx];}
        _sp_hybrid const int& dir() const {return (*this)[dir_idx];}
        _sp_hybrid const int& i  () const {return (*this)[i_idx];}
        _sp_hybrid const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        _sp_hybrid const int& j  () const {return (*this)[j_idx];}
        _sp_hybrid const int& k  () const {return (*this)[k_idx];}
        _sp_hybrid const int& lb () const {return (*this)[lb_idx];}
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
        
        _sp_hybrid node_idx_t(){}
        _sp_hybrid int& i () {return (*this)[i_idx];}
        _sp_hybrid int& i (const int& idx) {return (*this)[i_idx + idx];}
        _sp_hybrid int& j () {return (*this)[j_idx];}
        _sp_hybrid int& k () {return (*this)[k_idx];}
        _sp_hybrid int& lb() {return (*this)[lb_idx];}
        _sp_hybrid const int& i () const {return (*this)[i_idx];}
        _sp_hybrid const int& i (const int& idx) const {return (*this)[i_idx + idx];}
        _sp_hybrid const int& j () const {return (*this)[j_idx];}
        _sp_hybrid const int& k () const {return (*this)[k_idx];}
        _sp_hybrid const int& lb() const {return (*this)[lb_idx];}
    };
    
    //i_cell is the cell index
    //idir is the principal cartesian component of the face normal vector [unsigned]
    //pm = 0 indicates the negative-side face, pm = 1 indicates the positive-side face
    _sp_hybrid static face_idx_t cell_to_face(const cell_idx_t& i_cell, const int& idir, const int& pm)
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
    _sp_hybrid static node_idx_t cell_to_node(const cell_idx_t& i_cell, const ctrs::array<int,3>& offset)
    {
        node_idx_t output;
        output.i()  = i_cell.i();
        output.j()  = i_cell.j();
        output.k()  = i_cell.k();
        output.lb() = i_cell.lb();
        output.i()  += offset[0];
        output.j()  += offset[1];
        output.k()  += offset[2];
        return output;
    }
    
    //i_face is the face index
    //pm = 0 indicates the negative-side ("left") face, pm = 1 indicates the positive-side ("right") face
    _sp_hybrid static cell_idx_t face_to_cell(const face_idx_t& i_face, const int& pm)
    {
        cell_idx_t output;
        output.i() =  i_face.i();
        output.j() =  i_face.j();
        output.k() =  i_face.k();
        output.lb() = i_face.lb();
        output.i(i_face.dir()) += (pm-1);

        return output;
    }
    
    //i_face is the face index
    //pm = 0 indicates the negative-side ("left") face, pm = 1 indicates the positive-side ("right") face
    _sp_hybrid static node_idx_t face_to_node(const face_idx_t& i_face, const int& pm)
    {
        // print("NOT IMPLEMENTED:", __FILE__, __LINE__);
        abort();
        node_idx_t output((int)i_face[0], (int)i_face[1], (int)i_face[2], (int)i_face[3]);
        return output;
    }
    
    //Return the index-space coordinates of the index, used to compute the computational coordinates
    template <typename rtype=real_t> _sp_hybrid ctrs::array<rtype, 3> get_index_coord(const cell_idx_t& i_cell)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_cell.i()+0.5,
            (rtype)i_cell.j()+0.5,
            (rtype)i_cell.k()+0.5);
        return output;
    }
    
    template <typename rtype=real_t> _sp_hybrid ctrs::array<rtype, 3> get_index_coord(const node_idx_t& i_node)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_node.i(),
            (rtype)i_node.j(),
            (rtype)i_node.k());
        return output;
    }
    
    template <typename rtype=real_t> _sp_hybrid ctrs::array<rtype, 3> get_index_coord(const face_idx_t& i_face)
    {
        ctrs::array<rtype, 3> output(
            (rtype)i_face.i()+0.5,
            (rtype)i_face.j()+0.5,
            (rtype)i_face.k()+0.5);
        output[i_face.dir()] -= 0.5;
        return output;
    }
    
    template <const array_centering centering> struct get_index_type{};
    template <> struct get_index_type<cell_centered>
    {
        typedef typename cell_idx_t::value_type integral_type;
        typedef cell_idx_t array_type;
    };
    template <> struct get_index_type<face_centered>
    {
        typedef typename face_idx_t::value_type integral_type;
        typedef face_idx_t array_type;
    };
    template <> struct get_index_type<node_centered>
    {
        typedef typename node_idx_t::value_type integral_type;
        typedef node_idx_t array_type;
    };
    template <> struct get_index_type<edge_centered>
    {
        typedef typename node_idx_t::value_type integral_type;
        typedef edge_idx_t array_type;
    };
}