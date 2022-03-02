#pragma once
#include <string>
#include <concepts>

#include "ctrs.h"
namespace cvdf::fluid_state
{
    
    template <class T> concept is_state_type = requires(T t, size_t idx)
    {
        t[idx];
        T::size();
        t.name(0);
    };
    
    template <typename rtype> struct prim_t
    {
        static constexpr size_t size(void) noexcept {return 5;}
        rtype data[size()];
        rtype& p() {return data[0];}
        rtype& T() {return data[1];}
        rtype& u() {return data[2];}
        rtype& v() {return data[3];}
        rtype& w() {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
        std::string name(uint idx) const
        {
            ctrs::array<std::string, 5> names("P", "T", "U", "V", "W");
            return names[idx];
        }
    };

    template <typename rtype> struct cons_t
    {
        static constexpr size_t size(void) noexcept {return 5;}
        rtype data[size()];
        rtype& rho  () {return data[0];}
        rtype& rho_H() {return data[1];}
        rtype& rho_u() {return data[2];}
        rtype& rho_v() {return data[3];}
        rtype& rho_w() {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
        std::string name(uint idx) const
        {
            ctrs::array<std::string, 5> names("rho", "rhoH", "rhoU", "rhoV", "rhoW");
            return names[idx];
        }
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

    template<typename ptype, typename ctype, class gastype> void convert_state(prim_t<ptype>& prim, cons_t<ctype>& cons, const gastype& gas)
    {
        ptype rho = prim.p() / (gas.R*prim.T());
        ptype rhoU2 = rho*(prim.u()*prim.u()+prim.v()*prim.v()+prim.w()*prim.w());
        ptype rhoE = 0.5*rhoU2 + (prim.p()/((gas.gamma - 1.0)));
        ptype rhoU = rho*prim.u();
        ptype rhoV = rho*prim.v();
        ptype rhoW = rho*prim.w();
        cons.rho()   = rho;
        cons.rho_H() = rhoE;
        cons.rho_u()  = rhoU;
        cons.rho_v()  = rhoV;
        cons.rho_w()  = rhoW;
    }

    template<typename ptype, typename ctype, class gastype> void convert_state(cons_t<ctype>& cons, prim_t<ptype>& prim, const gastype& gas)
    {
        ptype rho = cons.rho();
        ptype invrho = 1.0/rho;
        ptype u = invrho*cons.rho_u();
        ptype v = invrho*cons.rho_v();
        ptype w = invrho*cons.rho_w();
        ptype rhoU2 = rho*(u*u+v*v+w*w);
        ptype p = (gas.gamma - 1.0)*(cons.rho_H() - 0.5*rhoU2);
        ptype T = p/(gas.R*rho);
        prim.p() = p;
        prim.T() = T;
        prim.u() = u;
        prim.v() = v;
        prim.w() = w;
    }
}