#pragma once

#include "core/cuda_incl.h"
#include "omni/omni.h"
#include "core/consts.h"
#include <fstream>
#include <string>
#include <vector>

namespace spade::fluid_state
{
	template <class T> concept is_state_type = requires(T t, size_t idx)
	{
		t[idx];
		T::size();
		t.name(0);
		typename T::value_type;
	};

	// This could probably be implemented better <JRB>
	template <class T> concept is_multicomponent_gas_type = requires(T t)
	{
		t.nspecies();
		t.mw_s[0]; t.mw_si[0]; t.hf_s[0]; t.charge_s[0];
		t.isMol[0];
	};

	template <typename derived_t, typename ilist_t = omni::info_list_t<>>
	struct gas_interface_t
	{
		_sp_hybrid derived_t&       self()       {return *static_cast<      derived_t*>(this);}
		_sp_hybrid const derived_t& self() const {return *static_cast<const derived_t*>(this);}

		using info_type = ilist_t;

        _sp_hybrid auto get_R    (const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_R    (args...);}, input);
        }

        _sp_hybrid auto get_gamma(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_gamma(args...);}, input);
        }
    };

	template <typename dtype> struct ideal_gas_t
	: public gas_interface_t<ideal_gas_t<dtype>>
	{
		using base_t = gas_interface_t<ideal_gas_t<dtype>>;
		using base_t::get_R;
		using base_t::get_gamma;
		using base_t::info_type;

		typedef dtype value_type;
		dtype R, gamma;
		_sp_hybrid ideal_gas_t(){}
		_sp_hybrid ideal_gas_t(const dtype& gamma_in, const dtype& R_in) : gamma{gamma_in}, R{R_in} {}
		_sp_hybrid dtype get_R    () const {return this->R;}
		_sp_hybrid dtype get_gamma() const {return this->gamma;}
	};

	template <typename dtype, const std::size_t num_species, const std::size_t maxVLevel> 
	struct multicomponent_gas_t : public gas_interface_t<multicomponent_gas_t<dtype, num_species, maxVLevel>>
	{
		using float_t = dtype;
		
		// Some member variables
		spade::ctrs::array<dtype, num_species> mw_s; // Molecular weight
		spade::ctrs::array<dtype, num_species> mw_si; // Inverse of molecular weight
		spade::ctrs::array<dtype, num_species> hf_s; // Heat of formation
		spade::ctrs::array<int, num_species> nvib; // Number of vibrationally activated energy levels
		spade::linear_algebra::dense_mat<dtype, num_species, maxVLevel> gvib; // Vibrational energy level degeneracy
		spade::linear_algebra::dense_mat<dtype, num_species, maxVLevel> theta_v; // Vibrational energy level characteristic temp.
		spade::ctrs::array<int, num_species>   isMol; // Molecule identification flag (0) atom, (1) molecule
		
		// Constructors
		_sp_hybrid multicomponent_gas_t(){}

		// Get species count
		_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}

		// Get max vibrational energy levels
		_sp_hybrid constexpr static std::size_t maxVibLevel(){return maxVLevel;}

		// Function -- compute species gas constant
		_sp_hybrid dtype get_Rs(const int s) const {return spade::consts::Rgas_uni * mw_si[s];}
		
		// Function -- compute translational specific heat
		_sp_hybrid dtype get_cvt(const int s) const {return float_t(1.5) * get_Rs(s);}

		// Function -- compute rotational specific heat
		_sp_hybrid dtype get_cvr(const int s) const {return get_Rs(s) * float_t(isMol[s]);}

		// Function -- compute translational/rotational specific heat
		_sp_hybrid dtype get_cvtr(const int s) const {return get_cvt(s) + get_cvr(s);}

		// Function -- compute vibrational specific heat
		_sp_hybrid dtype get_cvv(const int s, const dtype T) const
		{
			if (isMol[s]>0)
			{
				dtype Tinv = float_t(1.0) / T;
				dtype cvv = dtype(0.0);
				for (int m = 0; m<nvib[s]; ++m)
				{
					cvv += gvib(s,m) * get_Rs(s) * (theta_v(s,m) * Tinv) * (theta_v(s,m) * Tinv) * exp(theta_v(s,m) * Tinv) / ((exp(theta_v(s,m) * Tinv) - float_t(1.0)) * (exp(theta_v(s,m) * Tinv) - float_t(1.0)));
				}
				return cvv;
			}
			else
			{
				return 0.0;
			}
		}

		/*_sp_hybrid auto get_gamma(const auto& input) const
		{
			return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_gamma(args...);}, input);
		}*/
	};

	// Initialization function for incoming species vibrational energy data
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	static void import_vibrational_data(const std::string& vib_fname, const std::vector<std::string>& speciesNames, multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		std::ifstream infile;

		// Get spade environment path
		const auto env_p = std::getenv("SPADE");
		
		// Set full filename for vibrational data
		std::string full_fname = "";
		full_fname += env_p;
		full_fname += "/src/navier-stokes/speciesInputs/" + vib_fname;
		
		// Open file
		infile.open(full_fname);
		
		if (infile)
		{
			// String to read species name
			std::string species;

			// Temporary variables
			ptype gvib;
			ptype theta_v;
			gas.nvib = 0;
			
			// Sweep entire file
			int count = 0;
			while (true)
			{
				// Read line in file
				infile >> species >> gvib >> theta_v;
				
				// Sweep species
				for (int s = 0; s<gas.nspecies(); ++s)
				{
					// Do we need this species?
					if (species == speciesNames[s])
					{
						// Store data
						gas.gvib(s,gas.nvib[s])    = gvib;
						gas.theta_v(s,gas.nvib[s]) = theta_v;

						// Count energy levels
						++gas.nvib[s];
					}
				}
				// Check for end of file
				if (infile.eof()) break;
			}
			
			// Close file
			infile.close();
			
		}
		else
		{
			std::cerr << "Can not open provided species vibrational energy data file!" << std::endl;
		}
		
		return;
	}
	
	// Initialization function for incoming species data
	template<typename ptype, const std::size_t ns, const std::size_t nvib>
	static void import_species_data(const std::string& species_fname, const std::string& vib_fname, const std::vector<std::string>& speciesNames, multicomponent_gas_t<ptype, ns, nvib>& gas)
	{
		std::ifstream infile;

		// Get spade environment path
		const auto env_p = std::getenv("SPADE");

		// Set full filename
		std::string full_fname = "";
		full_fname += env_p;
		full_fname += "/src/navier-stokes/speciesInputs/" + species_fname;
		
		// Open file
		infile.open(full_fname);
		
		if (infile)
		{
			// String to read species name
			std::string species;

			// Temporary variables
			ptype mw, hf, charge;
			int isMol;
			
			// Sweep entire file
			int count = 0;
			while (true)
			{
				// Read line in file
				infile >> species >> mw >> isMol >> hf >> charge;
				
				// Sweep species
				for (int s = 0; s<gas.nspecies(); ++s)
				{
					// Do we need this species?
					if (species == speciesNames[s])
					{
						// Store data
						gas.mw_s[s]    = mw;
						gas.mw_si[s]   = float_t(1.0) / mw;
						gas.isMol[s]   = isMol;
						gas.hf_s[s]    = hf;
						gas.charge_s[s]  = charge;

						// Count species
						count += 1;
					}
				}
				// Check for end of file
				if (infile.eof()) break;
			}
			
			// Close file
			infile.close();
			
		}
		else
		{
			std::cerr << "Can not open provided species data file!" << std::endl;
		}

		// Run import on vibrational energy data now
		import_vibrational_data(vib_fname, speciesNames, gas);
		
		return;
	}
	
}
