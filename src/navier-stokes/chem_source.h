#pragma once
#include <string>

#include "core/config.h"
#include "core/ctrs.h"

#include "navier-stokes/gas.h"
#include "core/consts.h"
#include "navier-stokes/fluid_state.h"

namespace spade::fluid_state
{
	
	// Stores anything related to the employed reaction mechanism and/or vibrational relaxation model
	template<typename rtype> struct reactionMechanism_t
	{
		// Constructor
		_sp_hybrid reactionMechanism_t(){}

		// Number of reactions in mechanism
		static const int nr=5; // Hard-coded for 5-species air model for now

		// Reaction parameters
		constexpr static rtype Tmin=800;
		constexpr static rtype eps_chem2=6400; // 80^2
		spade::ctrs::array<int, 5> nIntervals; // Over species
		std::vector<rtype> Tlower,Tupper; // Dynamic vectors storing temp intervals for gibbs curve fit
		spade::ctrs::array<int, 5> reactType; // Over reactions
		spade::ctrs::array<rtype, 5> Af,Ta,eta; // Over reactions
		spade::ctrs::array<rtype, 9*5> gibbsCoefT1,gibbsCoefT2,gibbsCoefT3; // Over species (9 coef per species)
		spade::ctrs::array<rtype, 5*5> vprod,vreact; // Stoichiometric coefficients (ns by nr)
		spade::ctrs::array<rtype, 5> EII; // Electron impact ionization source term
		
		// Vibrational relaxation parameters
		constexpr static rtype p0=101325.0;
		constexpr static rtype Asr=0.00116,Bsr=0.015; // For vibrational relaxation		
		
	};
	
	// Overall chemical reaction mechanism structure
	template <typename rtype> struct chem_t
	: public spade::ctrs::arithmetic_array_t<rtype, 5, chem_t<rtype>>
	{
        using base_t = spade::ctrs::arithmetic_array_t<rtype, 5, chem_t<rtype>>;        
        using base_t::base_t;

		// Constructor
		_sp_hybrid chem_t(){}

		// Member variables
		_sp_hybrid rtype& omega(const int i) {return (*this)[i];}
		_sp_hybrid const rtype& omega(const int i) const {return (*this)[i];}

	};

	// Vibrational relaxation mechanism
	template <typename rtype> struct vib_t
	: public spade::ctrs::arithmetic_array_t<rtype, 4, vib_t<rtype>>
	{
		using base_t = spade::ctrs::arithmetic_array_t<rtype, 4, vib_t<rtype>>;
        using base_t::base_t;

		// Constructor
		_sp_hybrid vib_t(){}
		
		// Member variables
		_sp_hybrid rtype& St2v() {return (*this)[0];}
		_sp_hybrid rtype& Sc2v() {return (*this)[1];}
		_sp_hybrid rtype& Sh2e() {return (*this)[2];}
		_sp_hybrid rtype& Se2i() {return (*this)[3];}
		_sp_hybrid const rtype& St2v() const {return (*this)[0];}
		_sp_hybrid const rtype& Sc2v() const {return (*this)[1];}
		_sp_hybrid const rtype& Sh2e() const {return (*this)[2];}
		_sp_hybrid const rtype& Se2i() const {return (*this)[3];}
		
	};

	// Import gibbs energy data
	template<typename ptype>
	static void import_gibbsEnergy_data(const std::string& fname, const int& ns, const std::vector<std::string>& speciesNames, reactionMechanism_t<ptype>& react)
	{
		std::ifstream infile;
		try
		{
			infile.open(fname);
			
			// String to read species name
			std::string species;

			// Temporary variables
			spade::ctrs::array<int, 5> nIntervals = 0;
			ptype Tlower,Tupper;
			spade::ctrs::array<ptype, 9> gibbsCoef;
			
			// Sweep entire file
			while (true)
			{
				// Read line in file
				infile >> species >> Tlower >> Tupper >> gibbsCoef[0] >> gibbsCoef[1] >> gibbsCoef[2] >> gibbsCoef[3] >> gibbsCoef[4] >> gibbsCoef[5] >> gibbsCoef[6] >> gibbsCoef[7] >> gibbsCoef[8];
				
				// Sweep species
				for (int s = 0; s<ns; ++s)
				{
					// Do we need this species?
					if (species == speciesNames[s])
					{
						// Store lower/upper bound
						react.Tlower(nIntervals[s], s) = Tlower;
						react.Tupper(nIntervals[s], s) = Tupper;

						// Store coefficients
						if (nIntervals[s]==0)
						{
							for (int i = 0; i<gibbsCoef.size(); ++i) react.gibbsCoefT1(i,s) = gibbsCoef[i];
						}
						else if (nIntervals[s]==1)
						{
							for (int i = 0; i<gibbsCoef.size(); ++i) react.gibbsCoefT2(i,s) = gibbsCoef[i];
						}
						else if (nIntervals[s]==2)
						{
							for (int i = 0; i<gibbsCoef.size(); ++i) react.gibbsCoefT3(i,s) = gibbsCoef[i];
						}
						else
						{
							std::cerr << "Cannot hold more than 3 temperature intervals! Need to implement some 3D arrays maybe..." << std::endl;
						}
							
						// Count intervals
						nIntervals[s] += 1;
						print(species);
					}
				}
				// Check for end of file
				if (infile.eof()) break;
			}
			
			// Close file
			infile.close();
			
		}
		catch (...)
		{
			std::cerr << "Can not open provides reaction data file!" << std::endl;
		}

		return;
	}

	// Compute chemical source term
	template<typename ptype>
	_sp_hybrid static void compute_chemSource(chem_t<ptype>& source, const prim_chem_t<ptype>& prim, const multicomponent_gas_t<ptype>& gas, const reactionMechanism_t<ptype>& react)
	{
		// Allocate vectors for reactions

		// Compute forward reaction rates
		
		// Compute backward reaction rates
		

		// Compute reaction product
		for (int s = 0; s<prim.ns; ++s) source.omega(s) = 0.0;
		
		return;
	}

	// Compute translational to vibrational source term
	template<typename ptype>
	_sp_hybrid static void compute_St2v(vib_t<ptype>& sourceVib, const prim_chem_t<ptype>& prim, const multicomponent_gas_t<ptype>& gas, const reactionMechanism_t<ptype>& react)
	{
		// Compute density
		ptype rho = get_rho(prim, gas);
		spade::ctrs::array<ptype, prim.ns> tau_s = 1.0;

		// Vibrational energy
		spade::ctrs::array<ptype, prim.ns> ev_st = get_evs(prim.T(), gas);
		spade::ctrs::array<ptype, prim.ns> ev_s  = get_evs(prim.Tv(), gas);
			
		// Compute relaxation time

		// Compute source term component
		sourceVib.St2v() = 0.0;
		for (int s = 0; s<prim.ns; ++s)
		{
			if (gas.isMol[s]>0) sourceVib.St2v() += rho * prim.Ys(s) * (ev_st[s] - ev_s[s]) / tau_s[s];
		}
		return;
	}

	// Compute chemical to vibrational source term
	template<typename ptype>
    _sp_hybrid static void compute_Sc2v(vib_t<ptype>& sourceVib, const chem_t<ptype>& chemSource, const prim_chem_t<ptype>& prim, const multicomponent_gas_t<ptype>& gas, const reactionMechanism_t<ptype>& react)
	{
		// Vibrational energy
		spade::ctrs::array<ptype, prim.ns> ev_s  = get_evs(prim.Tv(), gas);
			
		// Compute source term component
		sourceVib.Sc2v() = 0.0;
		for (int s = 0; s<prim.ns; ++s)
		{
			if (gas.isMol[s]>0) sourceVib.Sc2v() += chemSource.omega(s) * ev_s[s];
		}
		return;
	}

	// Compute heavy particle to electron source term
	template<typename ptype>
    _sp_hybrid static void compute_Sh2e(vib_t<ptype>& sourceVib, const prim_chem_t<ptype>& prim, const multicomponent_gas_t<ptype>& gas)
	{
		// Compute source term component
		sourceVib.Sh2e() = 0.0;
		return;
	}

	// Compute electron impact ionization source term
	template<typename ptype>
	_sp_hybrid static void compute_Se2i(vib_t<ptype>& sourceVib, const chem_t<ptype>& chemSource, const prim_chem_t<ptype>& prim, const multicomponent_gas_t<ptype>& gas, const reactionMechanism_t<ptype>& react)
	{
		// Compute source term component
		sourceVib.Se2i() = 0.0;
		return;
	}

	// Chemical source term structure. Primary driver
	template<typename ptype> struct chem_source_t
	{
		using info_type   = omni::info_list_t<omni::info::value>;
		using omni_type   = omni::prefab::cell_mono_t<info_type>;
		using output_type = fluid_state::flux_chem_t<ptype>;
		using gas_t       = multicomponent_gas_t<ptype>;
		using react_t     = reactionMechanism_t<ptype>;

		// Stores multi-component gas model
		gas_t gas;

		// Stores reaction mechanism
		react_t reactions;

		chem_source_t(const gas_t& gas_in, const react_t& reactions_in) : gas{gas_in}, reactions{reactions_in} {}

		_sp_hybrid output_type operator() (const auto& input) const
		{
			// Access state vector
            const auto& q = omni::access<omni::info::value >(input.root());

			// Initialize out-going source term
			output_type output = 0.0;

			// Initialize structures
			chem_t<ptype> chemSource;
			vib_t<ptype> vibSource;

			// Compute chemical source term
			compute_chemSource(chemSource, q, gas, reactions);
			
			// Compute vibrational source term components
			compute_St2v(vibSource, q, gas, reactions);
			compute_Sc2v(vibSource, chemSource, q, gas, reactions);
			compute_Sh2e(vibSource, q, gas);
			compute_Se2i(vibSource, chemSource, q, gas, reactions);

			// Assign to out-going vector
			for (int s = 0; s<q.ns; ++s) output.continuity(s) = chemSource.omega(s);
			output.energyVib() = vibSource.St2v() + vibSource.Sc2v() + vibSource.Sh2e() - vibSource.Se2i();

			// Return vector for assignment onto the RHS
			return output;
		}
	};
	
}
