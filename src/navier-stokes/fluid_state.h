#pragma once
#include <string>

#include "core/config.h"
#include "core/ctrs.h"

#include "navier-stokes/gas.h"
#include "core/consts.h"

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

    template <typename rtype> struct prim_chem_t : public ctrs::arithmetic_array_t<rtype, 11, prim_chem_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 11, prim_chem_t<rtype>>;
        using base_t::base_t;
        static const int ns=5;
        _sp_hybrid prim_chem_t(){}
		_sp_hybrid rtype& YN2() {return (*this)[0];}
		_sp_hybrid rtype& YO2() {return (*this)[1];}
		_sp_hybrid rtype& YNO() {return (*this)[2];}
		_sp_hybrid rtype& YN() {return (*this)[3];}
		_sp_hybrid rtype& YO() {return (*this)[4];}
        _sp_hybrid rtype& Ys(const int i) {return (*this)[i];}
		_sp_hybrid rtype& p() {return (*this)[5];}
        _sp_hybrid rtype& T() {return (*this)[6];}
        _sp_hybrid rtype& Tv() {return (*this)[7];}
        _sp_hybrid rtype& u() {return (*this)[8];}
        _sp_hybrid rtype& u(const int i) {return (*this)[8+i];}
        _sp_hybrid rtype& v() {return (*this)[9];}
        _sp_hybrid rtype& w() {return (*this)[10];}
		_sp_hybrid const rtype& YN2() const {return (*this)[0];}
		_sp_hybrid const rtype& YO2() const {return (*this)[1];}
		_sp_hybrid const rtype& YNO() const {return (*this)[2];}
		_sp_hybrid const rtype& YN() const {return (*this)[3];}
		_sp_hybrid const rtype& YO() const {return (*this)[4];}
        _sp_hybrid const rtype& Ys(const int i) const {return (*this)[i];}
		_sp_hybrid const rtype& p() const {return (*this)[5];}
        _sp_hybrid const rtype& T() const {return (*this)[6];}
        _sp_hybrid const rtype& Tv() const {return (*this)[7];}
        _sp_hybrid const rtype& u() const {return (*this)[8];}
        _sp_hybrid const rtype& u(const int i) const {return (*this)[8+i];}
        _sp_hybrid const rtype& v() const {return (*this)[9];}
        _sp_hybrid const rtype& w() const {return (*this)[10];}
        static std::string name(uint idx)
        {
			ctrs::array<std::string, 11> names("YN2", "YO2", "YNO", "YN", "YO", "P", "T", "Tv", "U", "V", "W");
            return names[idx];
        }

    };

    template <typename rtype> struct cons_chem_t : public ctrs::arithmetic_array_t<rtype, 10, cons_chem_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 10, cons_chem_t<rtype>>;
        using base_t::base_t;
        static const int ns=5;
        _sp_hybrid cons_chem_t(){}
        _sp_hybrid rtype& rhos(const int i) {return (*this)[i];}
        _sp_hybrid rtype& E() {return (*this)[5];}
        _sp_hybrid rtype& Ev() {return (*this)[6];}
        _sp_hybrid rtype& rho_u() {return (*this)[7];}
        _sp_hybrid rtype& rho_v() {return (*this)[8];}
        _sp_hybrid rtype& rho_w() {return (*this)[9];}
		_sp_hybrid rtype& rho_u(const int i) {return (*this)[7+i];}
        _sp_hybrid const rtype& rhos(const int i) const {return (*this)[i];}
        _sp_hybrid const rtype& E() const {return (*this)[5];}
        _sp_hybrid const rtype& Ev() const {return (*this)[6];}
        _sp_hybrid const rtype& rho_u() const {return (*this)[7];}
        _sp_hybrid const rtype& rho_v() const {return (*this)[8];}
        _sp_hybrid const rtype& rho_w() const {return (*this)[9];}
		_sp_hybrid const rtype& rho_u(const int i) const {return (*this)[7+i];}
        static std::string name(uint idx)
        {
  	    ctrs::array<std::string, 10> names("rhoN2", "rhoO2", "rhoNO", "rhoN", "rhoO", "E", "Ev", "rhoU", "rhoV", "rhoW");
            return names[idx];
        }
      
    };
    
    template <typename rtype> struct flux_chem_t : public ctrs::arithmetic_array_t<rtype, 10, flux_chem_t<rtype>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 10, flux_chem_t<rtype>>;
        using base_t::base_t;
        static const int ns=5;
        _sp_hybrid flux_chem_t(){}
        _sp_hybrid rtype& continuity(const int i) {return (*this)[i];}
        _sp_hybrid rtype& energy    () {return (*this)[5];}
        _sp_hybrid rtype& energyVib () {return (*this)[6];}
        _sp_hybrid rtype& x_momentum() {return (*this)[7];}
        _sp_hybrid rtype& y_momentum() {return (*this)[8];}
        _sp_hybrid rtype& z_momentum() {return (*this)[9];}
        _sp_hybrid const rtype& continuity(const int i) const {return (*this)[i];}
        _sp_hybrid const rtype& energy    () const {return (*this)[5];}
        _sp_hybrid const rtype& energyVib () const {return (*this)[6];}
        _sp_hybrid const rtype& x_momentum() const {return (*this)[7];}
        _sp_hybrid const rtype& y_momentum() const {return (*this)[8];}
        _sp_hybrid const rtype& z_momentum() const {return (*this)[9];}
        static std::string name(uint idx)
        {
  	    ctrs::array<std::string, 10> names("continuityN2", "continuityO2", "continuityNO", "continuityN", "continuityO", "energy", "energyVib", "x_momentum", "y_momentum", "z_momentum");
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

	//
	// Some reacting flow wizardry
	//
	
	// Function -- compute mixture density (from primitive vector)
	template<typename ptype, class gas_t>
	_sp_hybrid static ptype get_rho(const prim_chem_t<ptype>& prim, const gas_t& gas)
	{
		ptype aux = 0.0;
		for (int s = 0; s<prim.ns; ++s)
		{
			if (gas.mw_s[s]>1)
			{
				aux += prim.Ys(s) * gas.get_Rs(s) * prim.T();
			}
			else
			{
				aux += prim.Ys(s) * gas.get_Rs(s) * prim.Tv();
			}
		}

		// Return density
		return prim.p() / aux;
	}

	// Function -- compute mixture density (from conservative vector)
	template<typename ptype>
	_sp_hybrid static ptype get_rho(const cons_chem_t<ptype>& cons)
	{
		ptype rho = 0.0;
		for (int s = 0; s<cons.ns; ++s) rho += cons.rhos(s);

		// Return density
		return rho;
	}
	
	// Function -- compute vibrational energy
	template<typename dtype, class gas_t>
	_sp_hybrid static spade::ctrs::array<dtype, 5> get_evs(const dtype& T, const gas_t& gas)
	{
		// Initialize to 0
		spade::ctrs::array<dtype, 5> ev_s = 0.0;
		for (int s=0; s<5; ++s)
		{
			// Check for molecules
			if (gas.isMol[s]>0) ev_s[s] += gas.get_Rs(s) * gas.theta_v[s] / (exp(gas.theta_v[s]/T) - float_t(1.0));
		}
		return ev_s;
	}

	// Function -- compute vibrational energy
	template<typename ptype, class gas_t>
	_sp_hybrid static ptype get_Ev(const prim_chem_t<ptype>& prim, const gas_t& gas)
	{
		// Initialize to 0
        ptype vibenergy = 0.0;
		ptype rho = get_rho(prim, gas);
		for (int s=0; s<prim.ns; ++s)
		{
			// Check for molecules
			if (gas.isMol[s]>0) vibenergy += rho * prim.Ys(s) * gas.get_Rs(s) * gas.theta_v[s] / (exp(gas.theta_v[s]/prim.Tv()) - float_t(1.0));
		}
		return vibenergy;
	}
	
	// Function -- compute total energy
	template<typename ptype, class gas_t>
	_sp_hybrid static ptype get_E(const prim_chem_t<ptype>& prim, const gas_t& gas)
	{
		// Compute kinetic energy first so we don't need to initialize to 0
		ptype rho = get_rho(prim, gas);
		ptype E   = float_t(0.5) * rho * (prim.u()*prim.u() + prim.v()*prim.v() + prim.w()*prim.w());
		spade::ctrs::array<ptype, prim.ns> ev_s=get_evs(prim.Tv(), gas);
		
		// Add internal energy contribution from each species
		for (int s=0; s<prim.ns; ++s) E += rho * prim.Ys(s) * (gas.get_cvtr(s) * prim.T() + gas.hf_s[s] + ev_s[s]);

		// Return value
		return E;
	}

	// Function -- compute multi-component speed of sound
	template<typename ptype, class gas_t>
	_sp_hybrid static ptype get_sos(const prim_chem_t<ptype>& prim, const gas_t& gas)
	{
		// Initialize
		ptype Cvtr = 0.0;

		// Compute sum(Y_s * Cvtr) = Cvtr
		for (int s = 0; s<prim.ns; ++s) Cvtr += prim.Ys(s) * gas.get_cvtr(s);

		// Compute beta
		ptype beta = 0.0;
		for (int s = 0; s<prim.ns; ++s) beta += prim.Ys(s) * gas.mw_si[s];
		beta *= spade::consts::Rgas_uni / Cvtr;

		// Compute speed of sound
		return sqrt((float_t(1.0) + beta) * prim.p() / spade::fluid_state::get_rho(prim, gas));
	};

	// Lets overload some functions. prim2cons transformation
    template<typename ptype, typename ctype, class gas_t>
    _sp_inline _sp_hybrid static void convert_state(const prim_chem_t<ptype>& prim, cons_chem_t<ctype>& cons, const gas_t& gas)
    {
		ptype rho     = get_rho(prim, gas);
		for (int s=0; s<prim.ns; ++s) cons.rhos(s)   = rho * prim.Ys(s);
		for (int i=0; i<3; ++i) cons.rho_u(i) = rho * prim.u(i);
		cons.E()      = get_E(prim, gas);
		cons.Ev()     = get_Ev(prim, gas);
    }

	// Lets overload some functions. cons2prim transformation
    template<typename ctype, typename ptype, class gas_t>
    _sp_inline _sp_hybrid static void convert_state(const cons_chem_t<ctype>& cons, prim_chem_t<ptype>& prim, const gas_t& gas)
    {
		ptype irho = float_t(1.0) / get_rho(cons);
		for (int s=0; s<prim.ns; ++s) prim.Ys(s)  = cons.rhos(s) * irho;
		for (int i=0; i<3; ++i) prim.u(i) = irho * cons.rho_u(i);

		// Some intermediate quantities (formation energy & translational/rotational energy)
		ptype hof = 0.0, rhoCv = 0.0;
		for (int s = 0; s<prim.ns; ++s)
		{
			hof   += cons.rhos(s) * gas.hf_s[s];
			rhoCv += cons.rhos(s) * gas.get_cvtr(s);
		}
		
		// Back out temperature
		prim.T()   = (cons.E() - float_t(0.5) * (cons.rho_u() * prim.u() + cons.rho_v() * prim.v() + cons.rho_w() * prim.w()) - cons.Ev() - hof) / rhoCv;

		// Get vibrational temeprature from non-linear Newton solve
		ptype Tv   = nonlinear_solve_Tv(prim, cons, gas);
		prim.Tv()  = Tv;
    }

	// Non-linear solver to get Tv from provided Ev
	template<typename ptype, typename ctype, class gas_t>
	_sp_hybrid static ptype nonlinear_solve_Tv(const prim_chem_t<ptype>& prim, const cons_chem_t<ctype>& cons, const gas_t& gas)
	{
		// Error tolerance on convergence
		ptype TOL=1E-5;
		int iterMax=100;

		// Initial guess on Tv
		prim_chem_t<ptype> prim_loc=prim;
		ptype Ev,dEv_dTv;
		prim_loc.Tv() = prim_loc.T();

		// Run iteration loop
		for (int iter=0; iter<iterMax; ++iter)
		{
			// Compute vibrational energy
			Ev      = get_Ev(prim_loc, gas);
			dEv_dTv = 0.0;
			for (int s = 0; s<prim.ns; ++s) dEv_dTv += cons.rhos(s) * gas.get_cvv(s,prim_loc.Tv());
			
			// Newton solver
			prim_loc.Tv() += (cons.Ev() - Ev) / dEv_dTv;
		}

		// Return computed vibrational temperature
		return prim_loc.Tv();
	}
	
}
