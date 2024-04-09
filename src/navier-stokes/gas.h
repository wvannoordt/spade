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


  // Reacting flow gas model
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

    template <typename dtype> struct multicomponent_gas_t
    {
	
		// Some member variables
		spade::ctrs::array<dtype, 5> mw_s,mw_si,hf_s,theta_v;
		spade::ctrs::array<int, 5>   isMol;
      
		// Constructors
		_sp_hybrid multicomponent_gas_t(){}

		// Function -- compute species gas constant
		_sp_hybrid dtype get_Rs(const int& s) const {return spade::consts::Rgas_uni * mw_si[s];}
      
		// Function -- compute translational specific heat
		_sp_hybrid dtype get_cvt(const int& s) const {return float_t(1.5) * get_Rs(s);}

		// Function -- compute rotational specific heat
		_sp_hybrid dtype get_cvr(const int& s) const {return get_Rs(s) * float_t(isMol[s]);}

		// Function -- compute translational/rotational specific heat
		_sp_hybrid dtype get_cvtr(const int& s) const {return get_cvt(s) + get_cvr(s);}

		// Function -- compute vibrational specific heat
		_sp_hybrid dtype get_cvv(const int& s, const dtype& T) const
		{
			if (isMol[s]>0)
			{
				dtype Tinv = float_t(1.0) / T;
				return get_Rs(s) * (theta_v[s] * Tinv) * (theta_v[s] * Tinv) * exp(theta_v[s] * Tinv) / ((exp(theta_v[s] * Tinv) - float_t(1.0)) * (exp(theta_v[s] * Tinv) - float_t(1.0)));
			}
			else
			{
				return 0.0;
			}
		}

    };

	// Initialization function for incoming species data
	template<typename ptype>
	static void import_species_data(const std::string& fname, const int& ns, const std::vector<std::string>& speciesNames, multicomponent_gas_t<ptype>& gas)
	{
		std::ifstream infile;
		try
		{
			infile.open(fname);

			// String to read species name
			std::string species;

			// Temporary variables
			ptype mw,hf,theta_v;
			int isMol;
			
			// Sweep entire file
			int count = 0;
			while (true)
			{
				// Read line in file
				infile >> species >> mw >> isMol >> hf >> theta_v;
				
				// Sweep species
				for (int s = 0; s<ns; ++s)
				{
					// Do we need this species?
					if (species == speciesNames[s])
					{
						// Store data
						gas.mw_s[s]    = mw;
						gas.mw_si[s]   = float_t(1.0) / mw;
						gas.isMol[s]   = isMol;
						gas.hf_s[s]    = hf;
						gas.theta_v[s] = theta_v;
						print(species);
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
		catch (...)
		{
			std::cerr << "Can not open provides species data file!" << std::endl;
		}

		return;
	}
	
}
