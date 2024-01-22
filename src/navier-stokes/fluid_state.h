#pragma once
#include <string>

#include "core/config.h"
#include "core/ctrs.h"

#include "navier-stokes/gas.h"

namespace spade::fluid_state
{   
    template <typename rtype> struct prim_t : public ctrs::arithmetic_array_t<rtype, 5, prim_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 5, prim_t<rtype>>;        
        using base_t::base_t;
        _sp_hybrid prim_t(){}
        _sp_hybrid rtype& p() {return (*this)[0];}
        _sp_hybrid rtype& T() {return (*this)[1];}
        _sp_hybrid rtype& u() {return (*this)[2];}
        _sp_hybrid rtype& u(const int i) {return (*this)[2+i];}
        _sp_hybrid rtype& v() {return (*this)[3];}
        _sp_hybrid rtype& w() {return (*this)[4];}
        _sp_hybrid const rtype& p() const {return (*this)[0];}
        _sp_hybrid const rtype& T() const {return (*this)[1];}
        _sp_hybrid const rtype& u() const {return (*this)[2];}
        _sp_hybrid const rtype& u(const int i) const {return (*this)[2+i];}
        _sp_hybrid const rtype& v() const {return (*this)[3];}
        _sp_hybrid const rtype& w() const {return (*this)[4];}
        static std::string name(uint idx)
        {
            ctrs::array<std::string, 5> names("P", "T", "U", "V", "W");
            return names[idx];
        }
    };

    template <typename rtype> struct cons_t : public ctrs::arithmetic_array_t<rtype, 5, cons_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 5, cons_t<rtype>>;
        using base_t::base_t;
        _sp_hybrid cons_t(){}
        _sp_hybrid rtype& rho  () {return (*this)[0];}
        _sp_hybrid rtype& rho_H() {return (*this)[1];}
        _sp_hybrid rtype& rho_u() {return (*this)[2];}
        _sp_hybrid rtype& rho_v() {return (*this)[3];}
        _sp_hybrid rtype& rho_w() {return (*this)[4];}
        _sp_hybrid const rtype& rho  () const {return (*this)[0];}
        _sp_hybrid const rtype& rho_H() const {return (*this)[1];}
        _sp_hybrid const rtype& rho_u() const {return (*this)[2];}
        _sp_hybrid const rtype& rho_v() const {return (*this)[3];}
        _sp_hybrid const rtype& rho_w() const {return (*this)[4];}
        static std::string name(uint idx)
        {
            ctrs::array<std::string, 5> names("rho", "rhoH", "rhoU", "rhoV", "rhoW");
            return names[idx];
        }
    };
    
    template <typename rtype> struct flux_t : public ctrs::arithmetic_array_t<rtype, 5, flux_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 5, flux_t<rtype>>;
        using base_t::base_t;
        _sp_hybrid flux_t(){}
        _sp_hybrid rtype& continuity() {return (*this)[0];}
        _sp_hybrid rtype& energy    () {return (*this)[1];}
        _sp_hybrid rtype& x_momentum() {return (*this)[2];}
        _sp_hybrid rtype& y_momentum() {return (*this)[3];}
        _sp_hybrid rtype& z_momentum() {return (*this)[4];}
        _sp_hybrid const rtype& continuity() const {return (*this)[0];}
        _sp_hybrid const rtype& energy    () const {return (*this)[1];}
        _sp_hybrid const rtype& x_momentum() const {return (*this)[2];}
        _sp_hybrid const rtype& y_momentum() const {return (*this)[3];}
        _sp_hybrid const rtype& z_momentum() const {return (*this)[4];}
        static std::string name(uint idx)
        {
            ctrs::array<std::string, 5> names("continuity", "energy", "x_momentum", "y_momentum", "z_momentum");
            return names[idx];
        }
    };
    
    template <class T> concept state_dependent_gas = std::floating_point<typename T::value_type> && requires(T t, prim_t<typename T::value_type> s)
    {
        t.get_R(s);
        t.get_gamma(s);
    };
    template <class T> concept state_independent_gas = requires(T t)
    {
        t.get_R();
        t.get_gamma();
    };
    
    template <is_state_type state_type> static std::ostream & operator<<(std::ostream & os, const state_type& state)
    {
       os << "{";
       for (size_t i = 0; i < state_type::size(); i++)
       {
           os << state_type::name(i) << ":" << state.data[i];
           if (i<state_type::size()-1) os << ", ";
       }
       os << "}";
       return os;
    }

    template<typename ptype, typename ctype, class gas_t>
    _sp_inline _sp_hybrid static void convert_state(const prim_t<ptype>& prim, cons_t<ctype>& cons, const gas_t& gas)
    {
        ptype rho = prim.p() / (gas.get_R()*prim.T());
        ptype rhoU2 = rho*(prim.u()*prim.u()+prim.v()*prim.v()+prim.w()*prim.w());
        ptype rhoE = ctype(0.5)*rhoU2 + (prim.p()/((gas.get_gamma() - ctype(1.0))));
        ptype rhoU = rho*prim.u();
        ptype rhoV = rho*prim.v();
        ptype rhoW = rho*prim.w();
        cons.rho()    = rho;
        cons.rho_H()  = rhoE;
        cons.rho_u()  = rhoU;
        cons.rho_v()  = rhoV;
        cons.rho_w()  = rhoW;
    }

    template<typename ptype, typename ctype, class gas_t>
    _sp_inline _sp_hybrid static void convert_state(const cons_t<ctype>& cons, prim_t<ptype>& prim, const gas_t& gas)
    {
        ptype rho = cons.rho();
        ptype invrho = ctype(1.0)/rho;
        ptype u = invrho*cons.rho_u();
        ptype v = invrho*cons.rho_v();
        ptype w = invrho*cons.rho_w();
        
        ptype rhoU2 = rho*(u*u+v*v+w*w);
        ptype p = (gas.get_gamma() - ctype(1.0))*(cons.rho_H() - ctype(0.5)*rhoU2);
        ptype T = p/(gas.get_R()*rho);
        prim.p() = p;
        prim.T() = T;
        prim.u() = u;
        prim.v() = v;
        prim.w() = w;
    }
    
    template<is_state_type stype, class gas_t>
    _sp_hybrid static void convert_state(const stype& in, stype& out, const gas_t& gas)
    {
        out = in;
    }
    
    template <typename from_t, typename to_t, typename gas_t> concept state_convertible = requires(const from_t& from, to_t& to, const gas_t& gas)
    {
        convert_state(from, to, gas);
    };
}