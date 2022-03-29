#pragma once
#include <string>
#include <concepts>

#include "core/ctrs.h"

namespace cvdf::fluid_state
{
    
    template <class T> concept is_state_type = requires(T t, size_t idx)
    {
        t[idx];
        T::size();
        t.name(0);
        typename T::value_type;
    };
    
    template <typename rtype> struct prim_t
    {
        typedef rtype value_type;
        static constexpr size_t size(void) noexcept {return 5;}
        rtype data[size()];
        prim_t(void){}
        prim_t(const rtype& vl){for (std::size_t i = 0; i < size(); i++) data[i]=vl;}
        rtype& p(void) {return data[0];}
        rtype& T(void) {return data[1];}
        rtype& u(void) {return data[2];}
        rtype& v(void) {return data[3];}
        rtype& w(void) {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
        
        const rtype& p(void) const {return data[0];}
        const rtype& T(void) const {return data[1];}
        const rtype& u(void) const {return data[2];}
        const rtype& v(void) const {return data[3];}
        const rtype& w(void) const {return data[4];}
        const rtype& operator [] (uint idx) const {return data[idx];}
        
        std::string name(uint idx) const
        {
            ctrs::array<std::string, 5> names("P", "T", "U", "V", "W");
            return names[idx];
        }
    };

    template <typename rtype> struct cons_t
    {
        typedef rtype value_type;
        static constexpr size_t size(void) noexcept {return 5;}
        rtype data[size()];
        cons_t(void){}
        cons_t(const rtype& vl){for (std::size_t i = 0; i < size(); i++) data[i]=vl;}
        rtype& rho  (void) {return data[0];}
        rtype& rho_H(void) {return data[1];}
        rtype& rho_u(void) {return data[2];}
        rtype& rho_v(void) {return data[3];}
        rtype& rho_w(void) {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
        
        const rtype& rho  (void) const {return data[0];}
        const rtype& rho_H(void) const {return data[1];}
        const rtype& rho_u(void) const {return data[2];}
        const rtype& rho_v(void) const {return data[3];}
        const rtype& rho_w(void) const {return data[4];}
        const rtype& operator [] (uint idx) const {return data[idx];}
        
        std::string name(uint idx) const
        {
            ctrs::array<std::string, 5> names("rho", "rhoH", "rhoU", "rhoV", "rhoW");
            return names[idx];
        }
    };
    
    template <class T> concept state_dependent_gas = requires(T t, prim_t<typename T::value_type> s)
    {
        t.get_R(s);
        t.get_gamma(s);
    };
    
    template <typename dtype> struct perfect_gas_t
    {
        typedef dtype value_type;
        dtype R, gamma;
        template <is_state_type state_t> dtype get_R(const state_t& state)     const {return this->R;}
        template <is_state_type state_t> dtype get_gamma(const state_t& state) const {return this->gamma;}
    };
    
    template <is_state_type state_type> static std::ostream & operator<<(std::ostream & os, const state_type& state)
    {
       os << "{";
       for (size_t i = 0; i < state_type::size(); i++)
       {
           os << state.name(i) << ":" << state.data[i];
           if (i<state_type::size()-1) os << ", ";
       }
       os << "}";
       return os;
    }

    template<typename ptype, typename ctype, class gas_t> void convert_state(prim_t<ptype>& prim, cons_t<ctype>& cons, const gas_t& gas)
    {
        ptype rho = prim.p() / (gas.get_R(prim)*prim.T());
        ptype rhoU2 = rho*(prim.u()*prim.u()+prim.v()*prim.v()+prim.w()*prim.w());
        ptype rhoE = 0.5*rhoU2 + (prim.p()/((gas.get_gamma(prim) - 1.0)));
        ptype rhoU = rho*prim.u();
        ptype rhoV = rho*prim.v();
        ptype rhoW = rho*prim.w();
        cons.rho()   = rho;
        cons.rho_H() = rhoE;
        cons.rho_u()  = rhoU;
        cons.rho_v()  = rhoV;
        cons.rho_w()  = rhoW;
    }

    template<typename ptype, typename ctype, class gas_t> void convert_state(cons_t<ctype>& cons, prim_t<ptype>& prim, const gas_t& gas)
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
}