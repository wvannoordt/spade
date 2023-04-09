#pragma once
#include <string>

#include "core/config.h"
#include "core/ctrs.h"

namespace spade::fluid_state
{
    
    template <class T> concept is_state_type = requires(T t, size_t idx)
    {
        t[idx];
        T::size();
        t.name(0);
        typename T::value_type;
    };
    
    template <typename rtype> struct prim_t : public ctrs::arithmetic_array_t<rtype, 5, prim_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 5, prim_t<rtype>>;        
        using base_t::base_t;
        prim_t(){}
        rtype& p() {return (*this)[0];}
        rtype& T() {return (*this)[1];}
        rtype& u() {return (*this)[2];}
        rtype& u(const int i) {return (*this)[2+i];}
        rtype& v() {return (*this)[3];}
        rtype& w() {return (*this)[4];}
        const rtype& p() const {return (*this)[0];}
        const rtype& T() const {return (*this)[1];}
        const rtype& u() const {return (*this)[2];}
        const rtype& u(const int i) const {return (*this)[2+i];}
        const rtype& v() const {return (*this)[3];}
        const rtype& w() const {return (*this)[4];}
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
        cons_t(){}
        rtype& rho  () {return (*this)[0];}
        rtype& rho_H() {return (*this)[1];}
        rtype& rho_u() {return (*this)[2];}
        rtype& rho_v() {return (*this)[3];}
        rtype& rho_w() {return (*this)[4];}
        const rtype& rho  () const {return (*this)[0];}
        const rtype& rho_H() const {return (*this)[1];}
        const rtype& rho_u() const {return (*this)[2];}
        const rtype& rho_v() const {return (*this)[3];}
        const rtype& rho_w() const {return (*this)[4];}
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
        flux_t(){}
        rtype& continuity() {return (*this)[0];}
        rtype& energy    () {return (*this)[1];}
        rtype& x_momentum() {return (*this)[2];}
        rtype& y_momentum() {return (*this)[3];}
        rtype& z_momentum() {return (*this)[4];}
        const rtype& continuity() const {return (*this)[0];}
        const rtype& energy    () const {return (*this)[1];}
        const rtype& x_momentum() const {return (*this)[2];}
        const rtype& y_momentum() const {return (*this)[3];}
        const rtype& z_momentum() const {return (*this)[4];}
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
    
    template <typename dtype> struct ideal_gas_t
    {
        typedef dtype value_type;
        dtype R, gamma;
        ideal_gas_t(){}
        ideal_gas_t(const dtype& gamma_in, const dtype& R_in) {gamma = gamma_in; R = R_in;}
        template <is_state_type state_t> dtype get_R(const state_t& state)     const {return this->R;}
        template <is_state_type state_t> dtype get_gamma(const state_t& state) const {return this->gamma;}
        dtype get_R()     const {return this->R;}
        dtype get_gamma() const {return this->gamma;}
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
    static void convert_state(const prim_t<ptype>& prim, cons_t<ctype>& cons, const gas_t& gas)
    {
        ptype rho = prim.p() / (gas.get_R(prim)*prim.T());
        ptype rhoU2 = rho*(prim.u()*prim.u()+prim.v()*prim.v()+prim.w()*prim.w());
        ptype rhoE = 0.5*rhoU2 + (prim.p()/((gas.get_gamma(prim) - 1.0)));
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
    static void convert_state(const cons_t<ctype>& cons, prim_t<ptype>& prim, const gas_t& gas)
    {
        ptype rho = cons.rho();
        ptype invrho = 1.0/rho;
        ptype u = invrho*cons.rho_u();
        ptype v = invrho*cons.rho_v();
        ptype w = invrho*cons.rho_w();
        ptype rhoU2 = rho*(u*u+v*v+w*w);
        ptype p = (gas.get_gamma(prim) - 1.0)*(cons.rho_H() - 0.5*rhoU2);
        ptype T = p/(gas.get_R(prim)*rho);
        prim.p() = p;
        prim.T() = T;
        prim.u() = u;
        prim.v() = v;
        prim.w() = w;
    }
    
    template<is_state_type stype, class gas_t> static void convert_state(const stype& in, stype& out, const gas_t& gas)
    {
        out = in;
    }
    
    template <typename from_t, typename to_t, typename gas_t> concept state_convertible = requires(const from_t& from, to_t& to, const gas_t& gas)
    {
        convert_state(from, to, gas);
    };
}