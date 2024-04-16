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
		_sp_hybrid rtype& rho_u(const int i) {return (*this)[2+i];}
		_sp_hybrid rtype& rho_v() {return (*this)[3];}
		_sp_hybrid rtype& rho_w() {return (*this)[4];}
		_sp_hybrid const rtype& rho  () const {return (*this)[0];}
		_sp_hybrid const rtype& rho_H() const {return (*this)[1];}
		_sp_hybrid const rtype& rho_u() const {return (*this)[2];}
		_sp_hybrid const rtype& rho_u(const int i) const {return (*this)[2+i];}
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

	template <typename rtype, const std::size_t num_species> struct prim_chem_t : public ctrs::arithmetic_array_t<rtype, 5+num_species, prim_chem_t<rtype, num_species>>
	{
		using base_t = ctrs::arithmetic_array_t<rtype, 5+num_species, prim_chem_t<rtype, num_species>>;
		using base_t::base_t;
	_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}
		_sp_hybrid prim_chem_t(){}
		_sp_hybrid rtype& Ys(const int i) {return (*this)[i];}
		_sp_hybrid rtype& p()  {return (*this)[num_species-1];}
		_sp_hybrid rtype& T()  {return (*this)[num_species];}
		_sp_hybrid rtype& Tv() {return (*this)[num_species+1];}
		_sp_hybrid rtype& u()  {return (*this)[num_species+2];}
		_sp_hybrid rtype& u(const int i)  {return (*this)[num_species+2+i];}
		_sp_hybrid rtype& v()  {return (*this)[num_species+3];}
		_sp_hybrid rtype& w()  {return (*this)[num_species+4];}
		_sp_hybrid const rtype& Ys(const int i) const {return (*this)[i];}
		_sp_hybrid const rtype& p() const  {return (*this)[num_species-1];}
		_sp_hybrid const rtype& T() const  {return (*this)[num_species];}
		_sp_hybrid const rtype& Tv() const {return (*this)[num_species+1];}
		_sp_hybrid const rtype& u() const  {return (*this)[num_species+2];}
		_sp_hybrid const rtype& u(const int i) const {return (*this)[num_species+2+i];}
		_sp_hybrid const rtype& v() const  {return (*this)[num_species+3];}
		_sp_hybrid const rtype& w() const  {return (*this)[num_species+4];}
		static std::string name(uint idx)
		{
		ctrs::array<std::string, 5+num_species> names;
		for (int n = 0; n<num_species-1; ++n) names[n] = "Y" + std::to_string(n);
		names[num_species-1] = "P";
		names[num_species  ] = "T";
		names[num_species+1] = "Tv";
		names[num_species+2] = "U";
		names[num_species+3] = "V";
		names[num_species+4] = "W";
			return names[idx];
		}

	};

	template <typename rtype, const std::size_t num_species> struct cons_chem_t : public ctrs::arithmetic_array_t<rtype, 5+num_species, cons_chem_t<rtype, num_species>>
	{
		using base_t = ctrs::arithmetic_array_t<rtype, 5+num_species, cons_chem_t<rtype, num_species>>;
		using base_t::base_t;
		_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}
        _sp_hybrid prim_chem_t(){}
        _sp_hybrid rtype& Ys(const int i) {return (*this)[i];}
		_sp_hybrid rtype& p() {return (*this)[num_species-1];}
		_sp_hybrid rtype& u() {return (*this)[num_species];}
        _sp_hybrid rtype& u(const int i) {return (*this)[num_species+i];}
        _sp_hybrid rtype& v() {return (*this)[num_species+1];}
        _sp_hybrid rtype& w() {return (*this)[num_species+2];}
        _sp_hybrid rtype& T() {return (*this)[num_species+3];}
        _sp_hybrid rtype& Tv() {return (*this)[num_species+4];}
        _sp_hybrid const rtype& Ys(const int i) const {return (*this)[i];}
		_sp_hybrid const rtype& p() const {return (*this)[num_species-1];}
        _sp_hybrid const rtype& u() const {return (*this)[num_species];}
        _sp_hybrid const rtype& u(const int i) const {return (*this)[num_species+i];}
        _sp_hybrid const rtype& v() const {return (*this)[num_species+1];}
        _sp_hybrid const rtype& w() const {return (*this)[num_species+2];}
		_sp_hybrid const rtype& T() const {return (*this)[num_species+3];}
        _sp_hybrid const rtype& Tv() const {return (*this)[num_species+4];}
        static std::string name(uint idx)
        {
			ctrs::array<std::string, 5+num_species> names;
			for (int n = 0; n<num_species-1; ++n) names[n] = "Y" + std::to_string(n);
			names[num_species-1] = "P";
			names[num_species  ] = "U";
			names[num_species+1] = "V";
			names[num_species+2] = "W";
			names[num_species+3] = "T";
			names[num_species+4] = "Tv";
            return names[idx];
        }

    };

    template <typename rtype, const std::size_t num_species> struct cons_chem_t : public ctrs::arithmetic_array_t<rtype, 5+num_species, cons_chem_t<rtype, num_species>>
    {
        using base_t = ctrs::arithmetic_array_t<rtype, 5+num_species, cons_chem_t<rtype, num_species>>;
        using base_t::base_t;
		_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}
        _sp_hybrid cons_chem_t(){}
        _sp_hybrid rtype& rhos(const int i) {return (*this)[i];}
		_sp_hybrid rtype& rho_u() {return (*this)[num_species];}
        _sp_hybrid rtype& rho_v() {return (*this)[num_species+1];}
        _sp_hybrid rtype& rho_w() {return (*this)[num_species+2];}
		_sp_hybrid rtype& rho_u(const int i) {return (*this)[num_species+i];}
        _sp_hybrid rtype& E() {return (*this)[num_species+3];}
        _sp_hybrid rtype& Ev() {return (*this)[num_species+4];}
        _sp_hybrid const rtype& rhos(const int i) const {return (*this)[i];}
		_sp_hybrid const rtype& rho_u() const {return (*this)[num_species];}
        _sp_hybrid const rtype& rho_v() const {return (*this)[num_species+1];}
        _sp_hybrid const rtype& rho_w() const {return (*this)[num_species+2];}
		_sp_hybrid const rtype& rho_u(const int i) const {return (*this)[num_species+i];}
        _sp_hybrid const rtype& E() const {return (*this)[num_species+3];}
        _sp_hybrid const rtype& Ev() const {return (*this)[num_species+4];}
        static std::string name(uint idx)
        {
			ctrs::array<std::string, 5+num_species> names;
			for (int n = 0; n<num_species; ++n) names[n] = "rho" + std::to_string(n);
			names[num_species  ] = "rhoU";
			names[num_species+1] = "rhoV";
			names[num_species+2] = "rhoW";
			names[num_species+3] = "E";
			names[num_species+4] = "Ev";
            return names[idx];
        }
      
    };
    
	template <typename rtype, const std::size_t num_species> struct flux_chem_t : public ctrs::arithmetic_array_t<rtype, 5+num_species, flux_chem_t<rtype, num_species>>
	{
		using base_t = ctrs::arithmetic_array_t<rtype, 5+num_species, flux_chem_t<rtype, num_species>>;
		using base_t::base_t;
		_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}
        _sp_hybrid flux_chem_t(){}
        _sp_hybrid rtype& continuity(const int i) {return (*this)[i];}
		_sp_hybrid rtype& x_momentum() {return (*this)[num_species];}
        _sp_hybrid rtype& y_momentum() {return (*this)[num_species+1];}
        _sp_hybrid rtype& z_momentum() {return (*this)[num_species+2];}
        _sp_hybrid rtype& energy    () {return (*this)[num_species+3];}
        _sp_hybrid rtype& energyVib () {return (*this)[num_species+4];}
        _sp_hybrid const rtype& continuity(const int i) const {return (*this)[i];}
		_sp_hybrid const rtype& x_momentum() const {return (*this)[num_species];}
        _sp_hybrid const rtype& y_momentum() const {return (*this)[num_species+1];}
        _sp_hybrid const rtype& z_momentum() const {return (*this)[num_species+2];}
        _sp_hybrid const rtype& energy    () const {return (*this)[num_species+3];}
        _sp_hybrid const rtype& energyVib () const {return (*this)[num_species+4];}
        static std::string name(uint idx)
        {
			ctrs::array<std::string, 5+num_species> names;
			for (int n = 0; n<num_species; ++n) names[n] = "continuity" + std::to_string(n);
			names[num_species] = "x-momentum";
			names[num_species+1] = "y-momentum";
			names[num_species+2] = "z-momentum";
			names[num_species+3] = "energy";
			names[num_species+4] = "energyVib";
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

	template <class T> concept is_flux_type = requires(T t)
	{
		t.get_R();
		t.get_gamma();
	};

    // Added to allow for generalization of different flux schemes and viscous models <JRB | Implemented 4-14-24>
	template <class T> concept is_prim_state_type = is_state_type<T> and requires(T t, size_t idx)
	{
        // ensures that we have access to following primitive variables
        t.p();
        t.T();
        t.u();
        t.v();
        t.w();
	};
	template <class T> concept is_cons_state_type = is_state_type<T> and requires(T t, size_t idx)
	{
        // ensures that we have access to following conservative variables
        t.rho();
        t.rho_H();
        t.rho_u();
        t.rho_v();
        t.rho_w();
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

	// Function -- compute species mass fractions
	template<typename ptype, const std::size_t ns>
	_sp_hybrid static spade::ctrs::array<ptype, ns> get_Ys(const prim_chem_t<ptype, ns>& prim)
	{
		// Initialize massfraction vector
		spade::ctrs::array<ptype, prim.nspecies()> Ys;

		// Copy 1->ns-1 massfractions over
		Ys[prim.nspecies()-1] = ptype(1.0);
		for (int s = 0; s<prim.nspecies()-1; s++)
		{
			Ys[s] = prim.Ys(s);
			Ys[prim.nspecies()-1] -= prim.Ys(s);
		}

		return Ys;
	}
	
	// Function -- compute species density
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> get_rhos(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{

		// Initialize vector
		spade::ctrs::array<ptype, prim.nspecies()> rhos;
		spade::ctrs::array<ptype, prim.nspecies()> Ys = get_Ys(prim);
		
		// Compute sum(Ys * Rs * T)
		ptype aux = ptype(0.0);
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			if (gas.mw_s[s]>1)
			{
				aux += Ys[s] * gas.get_Rs(s) * prim.T();
			}
			else
			{
				aux += Ys[s] * gas.get_Rs(s) * prim.Tv();
			}
		}

		// Compute density
		ptype rho = prim.p() / aux;

		// Compute species densities
		for (int s = 0; s<prim.nspecies(); ++s) rhos[s] = rho * Ys[s];

		return rhos;
	}

	// Function -- get pressure (for initial condition/BC)
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static ptype get_pressure(const ptype& rho, const spade::ctrs::array<ptype, ns>& Ys, const ptype& T, const ptype& Tv, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		ptype pressure = 0.0;
		for (int s = 0; s<Ys.size(); ++s)
		{
			if (gas.mw_s[s]>1)
			{
				pressure += rho * Ys[s] * gas.get_Rs(s) * T;
			}
			else
			{
				pressure += rho * Ys[s] * gas.get_Rs(s) * Tv;
			}
		}

		return pressure;
	}
	
	// Function -- molar fraction
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> get_Xr(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		using float_t = ptype;

		// Initialize
		spade::ctrs::array<ptype, prim.nspecies()> Xr;

		// Get species densities
		spade::ctrs::array<ptype, prim.nspecies()> rhos = get_rhos(prim, gas);
		
		// Compute molar fraction
		ptype sumXr = float_t(0.0);
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			sumXr += rhos[s] * gas.mw_si[s];
		}
		
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			Xr[s] = rhos[s] * gas.mw_si[s] / sumXr;
		}
		
		// Return output
		return Xr;
	}
	
	// Function -- compute mixture density (from primitive vector)
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static ptype get_rho(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		// Get species mass fractions
		spade::ctrs::array<ptype, prim.nspecies()> Ys = get_Ys(prim);
		

		// Compute sum(Ys * Rs * T)
		ptype aux = 0.0;
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			if (gas.mw_s[s]>1)
			{
				aux += Ys[s] * gas.get_Rs(s) * prim.T();
			}
			else
			{
				aux += Ys[s] * gas.get_Rs(s) * prim.Tv();
			}
		}

		// Return density
		return prim.p() / aux;
	}

	// Function -- compute mixture density (from conservative vector)
	template<typename ptype, const std::size_t ns>
	_sp_hybrid static ptype get_rho(const cons_chem_t<ptype, ns>& cons)
	{
		ptype rho = 0.0;
		for (int s = 0; s<cons.nspecies(); ++s) rho += cons.rhos(s);

		// Return density
		return rho;
	}
	
	// Function -- compute vibrational energy
	template<typename dtype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<dtype, ns> get_evs(const dtype& T, const multicomponent_gas_t<dtype, ns, nvib>& gas)
	{
		using float_t = dtype;
		
		// Initialize to 0
		spade::ctrs::array<dtype, gas.nspecies()> ev_s = 0.0;
		dtype Tinv = float_t(1.0) / T;
		for (int s=0; s<gas.nspecies(); ++s)
		{
			// Check for molecules
			if (gas.isMol[s]>0)
			{
				for (int m = 0; m<gas.nvib[s]; ++m)
				{
					ev_s[s] += gas.gvib(s,m) * gas.get_Rs(s) * gas.theta_v(s,m) / (exp(gas.theta_v(s,m)*Tinv) - float_t(1.0));
				}
			}
		}
		return ev_s;
	}

	// Function -- compute vibrational energy
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static ptype get_Ev(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		using float_t = ptype;

		// Get species densities
		spade::ctrs::array<ptype, prim.nspecies()> rhos = get_rhos(prim, gas);
		spade::ctrs::array<ptype, prim.nspecies()> ev_s = get_evs(prim.Tv(), gas);
		
		// Initialize to 0
        ptype vibenergy = 0.0;
		for (int s=0; s<prim.nspecies(); ++s)
		{
			// Check for molecules
			if (gas.isMol[s]>0) vibenergy += rhos[s] * ev_s[s];
		}
		
		return vibenergy;
	}
	
	// Function -- compute total energy
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static ptype get_E(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		using float_t = ptype;

		// Get species density and mixture density
		spade::ctrs::array<ptype, prim.nspecies()> rhos = get_rhos(prim, gas);
		ptype rho = 0.0;
		for (int s = 0; s<prim.nspecies(); ++s) rho += rhos[s];

		// Compute kinetic energy first so we don't need to initialize to 0
		ptype E   = float_t(0.5) * rho * (prim.u()*prim.u() + prim.v()*prim.v() + prim.w()*prim.w());
		spade::ctrs::array<ptype, prim.nspecies()> ev_s = get_evs(prim.Tv(), gas);
		
		// Add internal energy contribution from each species
		for (int s=0; s<prim.nspecies(); ++s) E += rhos[s] * (gas.get_cvtr(s) * prim.T() + gas.hf_s[s] + ev_s[s]);

		// Return value
		return E;
	}

	// Function -- compute multi-component speed of sound
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static ptype get_sos(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		using float_t = ptype;

		// Get species mass fraction
		spade::ctrs::array<ptype, prim.nspecies()> Ys = get_Ys(prim);
		
		// Initialize
		ptype Cvtr = 0.0;

		// Compute sum(Y_s * Cvtr_s) = Cvtr
		for (int s = 0; s<prim.nspecies(); ++s) Cvtr += Ys[s] * gas.get_cvtr(s);

		// Compute beta
		ptype beta = 0.0;
		for (int s = 0; s<prim.nspecies(); ++s) beta += Ys[s] * gas.mw_si[s];
		beta *= spade::consts::Rgas_uni / Cvtr;

		// Compute sum(Ys * Rs * T)
		ptype aux = 0.0;
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			if (gas.mw_s[s]>1)
			{
				aux += Ys[s] * gas.get_Rs(s) * prim.T();
			}
			else
			{
				aux += Ys[s] * gas.get_Rs(s) * prim.Tv();
			}
		}
		
		// Compute speed of sound
		return sqrt((float_t(1.0) + beta) * aux);
	};

	// Lets overload some functions. prim2cons transformation
    template<typename ptype, typename ctype, const std::size_t ns, const std::size_t nvib>
    _sp_inline _sp_hybrid static void convert_state(const prim_chem_t<ptype, ns>& prim, cons_chem_t<ctype, ns>& cons, const multicomponent_gas_t<ptype, ns, nvib>& gas)
    {
		// Get species density from primitive vector
		spade::ctrs::array<ptype, prim.nspecies()> rhos = get_rhos(prim, gas);

		// Compute mixture density
		ptype rho     = ptype(0.0);
		for (int s = 0; s<prim.nspecies(); ++s) rho += rhos[s];

		// Set conserved variables
		for (int s = 0; s<prim.nspecies(); ++s) cons.rhos(s)   = rhos[s];
		for (int i = 0; i<3; ++i) cons.rho_u(i) = rho * prim.u(i);
		cons.E()      = get_E(prim, gas);
		cons.Ev()     = get_Ev(prim, gas);
    }

	// Lets overload some functions. cons2prim transformation
    template<typename ctype, typename ptype, const std::size_t ns, const std::size_t nvib>
    _sp_inline _sp_hybrid static void convert_state(const cons_chem_t<ctype, ns>& cons, prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas)
    {
		using float_t = ptype;

		ptype irho = float_t(1.0) / get_rho(cons);
		for (int s=0; s<prim.nspecies(); ++s) prim.Ys(s)  = irho * utils::max(cons.rhos(s),float_t(1E-20));
		for (int i=0; i<3; ++i) prim.u(i) = irho * cons.rho_u(i);
		
		// Some intermediate quantities (formation energy & translational/rotational energy)
		ptype hof = float_t(0.0), rhoCv = float_t(0.0);
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			hof   += utils::max(cons.rhos(s),float_t(1E-20)) * gas.hf_s[s];
			rhoCv += utils::max(cons.rhos(s),float_t(1E-20)) * gas.get_cvtr(s);
		}
		
		// Back out temperature
		prim.T()   = (cons.E() - float_t(0.5) * (cons.rho_u() * prim.u() + cons.rho_v() * prim.v() + cons.rho_w() * prim.w()) - cons.Ev() - hof) / rhoCv;
		
		// Get vibrational temperature and pressure from non-linear Newton solve
		ptype Tv, p;
		nonlinear_solve(prim, cons, gas, Tv, p);
		prim.Tv()  = Tv;
		prim.p()   = p;
		
    }

	// Non-linear solver to get Tv from provided Ev
	template<typename ptype, typename ctype, const std::size_t ns, const std::size_t nvib>
	_sp_hybrid static void nonlinear_solve(const prim_chem_t<ptype, ns>& prim, const cons_chem_t<ctype, ns>& cons, const multicomponent_gas_t<ptype, ns, nvib>& gas, ptype& Tv, ptype& p)
	{
		// Error tolerance on convergence
		ptype TOL = ptype(1E-5);
		int iterMax = 100;

		// Initial guess on Tv
		prim_chem_t<ptype, prim.nspecies()> prim_loc=prim;
		ptype Ev,dEv_dTv;
		prim_loc.Tv() = prim_loc.T(); // Provide initial guess for Tv

		// Run iteration loop
		for (int iter=0; iter<iterMax; ++iter)
		{
			// Compute mixture pressure
			prim_loc.p() = ptype(0.0);
			for (int s = 0; s<prim.nspecies(); ++s)
			{
				if (gas.mw_s[s]>1)
				{
					prim_loc.p() += utils::max(cons.rhos(s),ptype(1E-20)) * gas.get_Rs(s) * prim_loc.T();
				}
				else
				{
					prim_loc.p() += utils::max(cons.rhos(s),ptype(1E-20)) * gas.get_Rs(s) * prim_loc.Tv();
				}
			}
			
			// Compute vibrational energy
			Ev      = get_Ev(prim_loc, gas);
			dEv_dTv = ptype(0.0);
			for (int s = 0; s<prim.nspecies(); ++s) dEv_dTv += spade::utils::max(cons.rhos(s), ptype(1E-20)) * gas.get_cvv(s,prim_loc.Tv());
			
			// Newton solver
			prim_loc.Tv() += (cons.Ev() - Ev) / dEv_dTv;

			// Tolerance check
			if (utils::abs(cons.Ev() - Ev) < TOL) break;
		}
		
		// Set returned values
		Tv = prim_loc.Tv();
		p  = prim_loc.p();
		
		return;
	}
	
}
