#pragma once
#include <string>

#include "core/config.h"
#include "core/ctrs.h"
#include "core/linear_algebra.h"

#include "navier-stokes/gas.h"
#include "core/consts.h"
#include "navier-stokes/fluid_state.h"

namespace spade::fluid_state
{
	
	// Reaction type enumeration
	enum reaction_type
	{
	 DISS_REACTION=1,
	 EXCH_REACTION=2,
	 EI_ION_REACTION=3
	};

	template<typename rtype, const std::size_t num_species, const std::size_t num_react, const bool is_image = false> struct reactionData_t
	{
		// Reaction parameters
		constexpr static rtype Tmin=800; // Minimum reaction evaluation temperature
		constexpr static rtype eps_chem2=6400; // 80^2 --> Used in reaction temperature limiting
		constexpr static rtype sigma_coef=1E-20; // Collisional cross-section for Park relaxation time
		constexpr static rtype kb_max=1E26; // Maximum reaction rate

		using image_type = reactionData_t<rtype, num_species, num_react, true>;

		template <typename data_t> using storage_t = typename std::conditional<is_image, utils::vec_image_t<data_t>, device::shared_vector<data_t>>::type;

		// Stores pointers to all reaction data
	    storage_t<int>   idata;
	    storage_t<rtype> rdata;

		// Integer data
		spade::utils::vec_image_t<int> nIntervals; // Over species
		spade::utils::vec_image_t<int> reactType; // Over reactions

		// Real data
		spade::utils::md_vec_image_t<rtype, 2, int> Tlower,Tupper; // Dynamic vectors storing temp intervals for gibbs curve fit
		spade::utils::vec_image_t<rtype> Af,Ta,eta; // Over reactions
		spade::utils::md_vec_image_t<rtype, 3, int> gibbsCoef; // Over species (9 coef per species)
		spade::utils::md_vec_image_t<rtype, 2, int> phi_diss; // Dissociation enhancement parameter
		spade::utils::md_vec_image_t<rtype, 2, int> vprod,vreact; // Stoichiometric coefficients
		
		// Vibrational relaxation parameters
		constexpr static rtype p0=101325.0;
		constexpr static rtype Asr=0.00116,Bsr=0.015; // For vibrational relaxation		
		
		reactionData_t()
		{
			if constexpr (!is_image)
			{
				// Count the amount of integer/real data in this structure
				const auto [isize, rsize] = assign_images(idata.data(device::cpu),rdata.data(device::cpu));

				// Re-size the pointer vectors to hold all the memory addresses
				idata.resize(isize);
				rdata.resize(rsize);

				// Re-compute pointer memory addresses as it likely changed during the re-size
				assign_images(idata.data(device::cpu),rdata.data(device::cpu));
			}
		}

		template <typename device_t>
		requires (!is_image)
		image_type image(const device_t& dev)
		{
			// Switching reactionData template to image type
			image_type output;

			// Make images of reaction data
			output.idata = utils::make_vec_image(idata.data(dev));
			output.rdata = utils::make_vec_image(rdata.data(dev));

			// Assign images
			output.assign_images(output.idata, output.rdata);
			return output;
		}

		void transfer()
		{
			// Transfer data onto GPU
			idata.transfer();
			rdata.transfer();
		}

		std::tuple<std::size_t, std::size_t> assign_images(auto& idata_in, auto& rdata_in)
		{
			// Initialize data counters
			std::size_t icount = 0;
			std::size_t rcount = 0;

			// Lambda to add data to storage vectors
			const auto add_vec = [&](auto& new_vec, const auto&... sizes)
			{
				// Is this a multi-dimensional vector?
				constexpr bool is_md = sizeof...(sizes) > 1;

				// Remove any "const" and "&" types from the data type to get base type
				using vec_t = typename utils::remove_all<decltype(new_vec)>::type;

				// Is this an integer type?
				constexpr bool is_integral = std::same_as<typename vec_t::value_type, int>;

				// Number of dimensions on multi-dimensional vector
				constexpr std::size_t rank = sizeof...(sizes);

				// Get total size of incoming vector. Product of all provided sizes
				new_vec.csize = (... * sizes);
				if constexpr (is_md) new_vec.sizes = {sizes...}; // Vector of all sizes

				if constexpr (is_integral)
				{
					// Integer counter and pointer assignment
					new_vec.ptr = &idata_in[0] + icount;
					icount += new_vec.csize;
				}
				else
				{
					// Real counter and pointer assignment
					new_vec.ptr = &rdata_in[0] + rcount;
					rcount += new_vec.csize;
				}
			};

			// Call lambda here and assign all parameters to structure
			add_vec(nIntervals, num_species);
			add_vec(reactType, num_react);
			add_vec(Tlower, num_species, 3);
			add_vec(Tupper, num_species, 3);
			add_vec(Af, num_react);
			add_vec(eta, num_react);
			add_vec(Ta, num_react);
			add_vec(gibbsCoef, num_species, 3, 9);
			add_vec(phi_diss, num_react, num_species);
			add_vec(vprod, num_react, num_species);
			add_vec(vreact, num_react, num_species);

			// Return total integer/real counters
			return std::make_tuple(icount, rcount);
		}
	};
	
	// Stores anything related to the employed reaction mechanism and/or vibrational relaxation model
	template<typename rtype, const std::size_t num_species, const std::size_t num_react, const bool is_image = false> struct reactionMechanism_t
	{
		using float_t = rtype;
		using image_type = reactionMechanism_t<rtype, num_species, num_react, true>;

		template <typename device_t>
		requires (!is_image)
		image_type image(const device_t& dev)
		{
			// Switching reactionData template to image type
			image_type output{reactionData.image(dev)};
			
			return output;
		}
		
		// Reaction parameters
		reactionData_t<rtype, num_species, num_react, is_image> reactionData;
		
		// Number of reactions in mechanism
		_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}
		_sp_hybrid constexpr static std::size_t nreact(){return num_react;}
		
		// Forward controlling temperature
		_sp_hybrid rtype compute_Tf(const int& r, const rtype& T, const rtype& Tv) const
		{
			// Check reaction type
			if (reactionData.reactType[r] == DISS_REACTION)
			{
				// Forward rate temperature for dissociation reactions averages T and Tv
				return sqrt(T * Tv); // Almost always follows this form. Worth it to generalize?
			}
			else
			{
				// Forward rate temperature always T for other reactions (unless we do ionization ...)
				return T;
			}
		}
		
		// Backward controlling temperature
		_sp_hybrid rtype compute_Tb(const int& r, const rtype& T, const rtype& Tv) const {return T;}
		
		// Reaction rate temperature limiting function
		_sp_hybrid rtype limiting(const rtype& T) const {return float_t(0.5) * ((T + reactionData.Tmin) + sqrt(pow(T-reactionData.Tmin,2) + reactionData.eps_chem2));}
		
		// Compute arrhenius curve fit
		_sp_hybrid rtype compute_rates(const int& r, const rtype& T) const {return reactionData.Af[r] * pow(T,reactionData.eta[r]) * exp(-reactionData.Ta[r] / T);}

		// Get necessary curve fitting coefficients for equilibrium constant
		_sp_hybrid void get_NASA9_coefficients(const int& s, const rtype& T, spade::ctrs::array<rtype, 9>& coefs) const
		{
			// Sweep temperature intervals
			for (int i = 0; i<reactionData.nIntervals[s]; ++i)
			{
				// Check intervals
				if (T >= reactionData.Tlower(s,i) && T < reactionData.Tupper(s,i))
				{
					// Select interval
					for (int n = 0; n<coefs.size(); ++n) coefs[n] = reactionData.gibbsCoef(s,i,n);
					return;
				}
			}

			// Handle extremities
			if (T <= reactionData.Tlower(s,0))
			{
				// First interval
				for (int n = 0; n<coefs.size(); ++n) coefs[n] = reactionData.gibbsCoef(s,0,n);
				return;
			}
			else if (T >= reactionData.Tupper(s,reactionData.nIntervals[s]-1))
			{
				// Last interval
				for (int n = 0; n<coefs.size(); ++n) coefs[n] = reactionData.gibbsCoef(s,reactionData.nIntervals[s]-1,n);
				return;
			}
			else
			{
				//std::cerr<<"How did you get here?? Invalid NASA9 coefficient access!"<<std::endl;
			}
		}

		// Evaluate enthalpy norm
		_sp_hybrid rtype evaluate_enthalpyNorm(const int& s, const rtype& T, const spade::ctrs::array<rtype, 9>& coefs) const
		{
			// Pre-compute some values
			rtype T2   = T*T;
			rtype T3   = T2*T;
			rtype T4   = T3*T;
			rtype logT = log(T);
			
			// Function evaluation
			return - coefs[0]/T2 + coefs[1] * logT/T + coefs[2] + coefs[3] * T * float_t(0.5) + coefs[4] * T2 * float_t(1.0)/float_t(3.0) + coefs[5] * T3 * float_t(0.25) + coefs[6] * T4 * float_t(0.2) + coefs[7]/T;
		}

		// Evaluate entropy norm
		_sp_hybrid rtype evaluate_entropyNorm(const int& s, const rtype& T, const spade::ctrs::array<rtype, 9>& coefs) const
		{
			// Pre-compute some values
			rtype T2   = T*T;
			rtype T3   = T2*T;
			rtype T4   = T3*T;
			rtype logT = log(T);

			// Function evaluation
			return -float_t(0.5) * coefs[0] / T2 - coefs[1] / T + coefs[2] * logT + coefs[3] * T + coefs[4] * T2 * float_t(0.5) + coefs[5] * T3 * float_t(1.0) / float_t(3.0) + coefs[6] * T4 * float_t(0.25) + coefs[8];
		}

		// Evaluate Landau-Teller inter-species relaxation time
		template<const std::size_t ns, const std::size_t nvib>
		_sp_hybrid rtype evaluate_tausr(const int& s, const int& r, const rtype& p, const rtype& T, const multicomponent_gas_t<rtype, ns, nvib>& gas) const
		{
			// Find theta min
			rtype theta_min = 1E9;
			for (int n = 0; n<gas.nvib[s]; ++n) theta_min = utils::min(theta_min,gas.theta_v(s,n));
			
			// Parameters
			rtype mu_sr = gas.mw_s[s] * gas.mw_s[r] / (gas.mw_s[s] + gas.mw_s[r]);
			rtype Acoef = reactionData.Asr * sqrt(mu_sr) * pow(theta_min,float_t(4.0)/float_t(3.0));
			rtype Bcoef = reactionData.Bsr * sqrt(sqrt(mu_sr));
			return (reactionData.p0 / p) * exp(Acoef * (pow(T,-float_t(1.0)/float_t(3.0)) - Bcoef) - float_t(18.42));
		}
		
	};

	// Import gibbs energy data
	template<typename ptype, const std::size_t ns, const std::size_t nr, const bool is_image>
	static void import_gibbsEnergy_data(const std::string& gibbs_fname, const std::vector<std::string>& speciesNames, reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		std::ifstream infile;

		// Get spade environment path
		const auto env_p = std::getenv("SPADE");

		// Set full filename
		std::string full_fname = "";
		full_fname += env_p;
		full_fname += "/src/navier-stokes/reactionMechanisms/" + gibbs_fname;

		// Open file
		infile.open(full_fname);
		
		if (infile)
		{
			
			// String to read species name
			std::string species;

			// Temporary variables
			spade::ctrs::array<int, react.nspecies()> nIntervals = 0;
			ptype Tlower,Tupper;
			spade::ctrs::array<ptype, 9> gibbsCoef;
			
			// Sweep entire file
			while (true)
			{
				// Read line in file
				infile >> species >> Tlower >> Tupper >> gibbsCoef[0] >> gibbsCoef[1] >> gibbsCoef[2] >> gibbsCoef[3] >> gibbsCoef[4] >> gibbsCoef[5] >> gibbsCoef[6] >> gibbsCoef[7] >> gibbsCoef[8];

				// Check for end of file
				if (infile.eof()) break;
				
				// Sweep species
				for (int s = 0; s<react.nspecies(); ++s)
				{
					// Do we need this species?
					if (species == speciesNames[s])
					{
						// Store lower/upper bound
						react.reactionData.Tlower(s, nIntervals[s]) = Tlower;
						react.reactionData.Tupper(s, nIntervals[s]) = Tupper;

						// Store coefficients
						for (int i = 0; i<gibbsCoef.size(); ++i) react.reactionData.gibbsCoef(s,nIntervals[s],i) = gibbsCoef[i];
							
						// Count intervals
						nIntervals[s] += 1;
						
					}
				}
			}
			
			// Close file
			infile.close();
			
			// Copy intervals to reaction structure
			for (int s = 0; s<react.nspecies(); ++s) react.reactionData.nIntervals[s] = nIntervals[s];
			
		}
		else
		{
			std::cerr << "Can not open provided gibbs energy data file!" << std::endl;
		}
		
		return;
	}

	// Import reaction mechanism data data
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	static void import_reaction_data(const std::string& react_fname, const std::string& gibbs_fname, const std::vector<std::string>& speciesNames, const multicomponent_gas_t<ptype, ns, nvib>& gas, reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		using float_t = ptype;
		std::ifstream infile;
		
		// Get spade environment path
		const auto env_p = std::getenv("SPADE");

		// Set full filename
		std::string full_fname = "";
		full_fname += env_p;
		full_fname += "/src/navier-stokes/reactionMechanisms/" + react_fname;

		// Open file
		infile.open(full_fname);
		
		if (infile)
		{

			
			// String to read species name
			std::string reactType,species;

			// Temporary variables
			int nEnhance;
			std::vector<std::string> participants(6);
			ptype Af,Ta,eta,phi;

			// Add reaction check
			int count;
			int numReact=0;
			
			// Sweep entire file
			while (true)
			{
				// Read line in file
				infile >> reactType >> participants[0] >> participants[1] >> participants[2] >> participants[3] >> participants[4] >> participants[5] >> Af >> eta >> Ta >> nEnhance;

				// Check for end of file
				if (infile.eof()) break;
				
				// Do we need this reaction
				count = 0;
				for (int p = 0; p<participants.size(); ++p)
				{
					if (participants[p] == "--" || participants[p] == "M")
					{
						count += 1;
					}
					else
					{
						for (int s = 0; s<react.nspecies(); ++s)
						{
							if (participants[p] == speciesNames[s]) count += 1;
						}
					}
				}

				// Do we need this reaction? Only if all participants are present in the mixture
				if (count == participants.size())
				{
					// Save data
					react.reactionData.Af[numReact]  = Af;
					react.reactionData.eta[numReact] = eta;
					react.reactionData.Ta[numReact]  = Ta;

					// Set reaction type
					if (reactType == "diss")
					{
						react.reactionData.reactType[numReact] = DISS_REACTION;
					}
					else if (reactType == "exch")
					{
						react.reactionData.reactType[numReact] = EXCH_REACTION;
					}
					else if (reactType == "ei_ion")
					{
						react.reactionData.reactType[numReact] = EI_ION_REACTION;
					}
					else
					{
						std::cerr << "Undefined reaction type read in!" << std::endl;						
					}

					// Set stoichiometric coefficients (reactants)
					for (int p = 0; p<participants.size()/2; ++p)
					{
						for (int s = 0; s<react.nspecies(); ++s)
						{
							if (participants[p] == speciesNames[s]) react.reactionData.vreact(numReact, s) += float_t(1.0);
						}
					}
					
					// Set stoichiometric coefficients (products)
					for (int p = participants.size()/2; p<participants.size(); ++p)
					{
						for (int s = 0; s<react.nspecies(); ++s)
						{
							if (participants[p] == speciesNames[s]) react.reactionData.vprod(numReact, s) += float_t(1.0);
						}
					}

					// Sweep enhancement factor lines
					for (int n = 0; n<nEnhance; ++n)
					{
						infile >> species >> phi;
						
						// Check blanket statements
						if (species == "mol")
						{
							// Enhance dissociation rates for molecular collisions
							for (int s = 0; s<react.nspecies(); ++s)
							{
								if (gas.isMol[s]>0) react.reactionData.phi_diss(numReact,s) = phi;
							}
						}
						else if (species == "atoms")
						{
							// Enhance dissociation rates for atomic collisions
							for (int s = 0; s<react.nspecies(); ++s)
							{
								if (!gas.isMol[s]>0) react.reactionData.phi_diss(numReact,s) = phi;
							}
						}
						else
						{
							// Enhance dissociation rates for specific collisions
							for (int s = 0; s<react.nspecies(); ++s)
							{
								if (species == speciesNames[s]) react.reactionData.phi_diss(numReact,s) = phi;
							}
						}
					}
					
					// Count number of reactions as storage index
					numReact += 1;
				}
				else
				{
					// Sweep enhancement factor lines to skip them
					for (int n = 0; n<nEnhance; ++n)
					{
						infile >> species >> phi;
					}
				}
				
			}
			
			// Close file
			infile.close();
		}
		else
		{
			std::cerr << "Can not open provided reaction data file!" << std::endl;
		}

		// Import gibbs energy file
		import_gibbsEnergy_data(gibbs_fname, speciesNames, react);
		
		// Transfer reaction data onto GPU
		react.reactionData.transfer();
		
		return;
	}

	// Compute chemical source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_chemSource(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		using float_t = ptype;
		
		// Allocate vectors for reactions
		ptype kfb;
		spade::ctrs::array<ptype, react.nreact()> kf,kb;
		spade::ctrs::array<ptype, prim.nspecies()> source = float_t(0.0);
		
		// Initialize
		ptype Tf,Tb,hi,si,gibbs,Kc,vr;
		spade::ctrs::array<ptype, 9> coefs;
		
		// Loop reactions and compute forward/backward reaction rates
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Forward controlling temperature
			Tf = react.compute_Tf(r, prim.T(), prim.Tv());
			
			// Forward-backward rate temperature
			Tb = react.compute_Tb(r, prim.T(), prim.Tv());

			// Limit temperature lower bound
			Tf = react.limiting(Tf);
			Tb = react.limiting(Tb);

			// Compute rates
			kf[r]  = react.compute_rates(r, Tf);
			kfb    = react.compute_rates(r, Tb);

			// Initialize gibbs energy
			gibbs = float_t(0.0);

			// Initialize stoichiometric summation
			vr = float_t(0.0);

			// Sweep species
			for (int s = 0; s<prim.nspecies(); ++s)
			{
				if (utils::abs(react.reactionData.vprod(r,s) - react.reactionData.vreact(r,s)) > 0.5)
				{
					// Get curve fitting coefficients
					react.get_NASA9_coefficients(s, Tb, coefs);
					
					// Enthalpy norm
					hi = react.evaluate_enthalpyNorm(s, Tb, coefs);
					
					// Entropy norm
					si = react.evaluate_entropyNorm(s, Tb, coefs);
					
					// Gibbs energy summation
					gibbs += (react.reactionData.vprod(r,s) - react.reactionData.vreact(r,s)) * (hi - si);
					
					// Stoichiometric summation
					vr += react.reactionData.vprod(r,s) - react.reactionData.vreact(r,s);
				}
			}

			// Compute equilibrium constant
			Kc = exp(-gibbs) * pow(float_t(0.1) / (spade::consts::Rgas_uni * Tb * 1E-3),vr);

			// Compute backward reaction rate
			const auto kb_max = react.reactionData.kb_max;
			kb[r] = utils::min(kfb / Kc, kb_max);
		}

		// Initialize some variables
		spade::ctrs::array<ptype, react.nreact()> con   = float_t(0.0);
		spade::ctrs::array<ptype, react.nreact()> Rf_kf = float_t(1000.0);
		spade::ctrs::array<ptype, react.nreact()> Rb_kb = float_t(1000.0);

		// Compute species density
		spade::ctrs::array<ptype, prim.nspecies()> rhos = fluid_state::get_rhos(prim, gas);

		// Sweep reactions
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Check for dissociation reactions
			if (react.reactionData.reactType[r] == DISS_REACTION)
			{
				// Sweep species
				for (int s = 0; s<prim.nspecies(); ++s)
				{
					// Block out electrons
					if (gas.mw_s[s]>1) con[r] += rhos[s] * react.reactionData.phi_diss(r,s) * gas.mw_si[s] * 1E-3;
				}
			}
			else
			{
				// Turn off dissociation enhancement
				con[r] = float_t(1.0);
			}
		}

		// Form reaction product
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Sweep species
			for (int s = 0; s<prim.nspecies(); ++s)
			{
				// Forward reaction rate
				if (react.reactionData.vreact(r,s) > float_t(0.5)) Rf_kf[r] *= pow(1E-3 * rhos[s] * gas.mw_si[s], react.reactionData.vreact(r,s));

				// Backward reaction rate
				if (react.reactionData.vprod(r,s) > float_t(0.5)) Rb_kb[r] *= pow(1E-3 * rhos[s] * gas.mw_si[s], react.reactionData.vprod(r,s));
			}
		}

		// Final computation for source term
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Sweep species in reaction
			for (int s = 0; s<prim.nspecies(); ++s)
			{
				source[s] += (react.reactionData.vprod(r,s) - react.reactionData.vreact(r,s)) * (kf[r] * Rf_kf[r] - kb[r] * Rb_kb[r]) * gas.mw_s[s] * con[r];
			}
		}
		
		return source;
	}

	// Computation for Park vibrational relaxation time for high temperature corrections
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_parkRelaxTime(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		using float_t = ptype;
		
		// Initialize relaxation time
		spade::ctrs::array<ptype, prim.nspecies()> tau = float_t(0.0);
		
		// Compute collision cross-section
		ptype sigma = react.reactionData.sigma_coef * float_t(2.5E9) / (prim.T() * prim.T());

		// Get species densities
		spade::ctrs::array<ptype, prim.nspecies()> rhos = fluid_state::get_rhos(prim, gas);

		// Loop molecules
		ptype cs, Num_den;
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			if (gas.isMol[s]>0)
			{
				// Molecular velocity
				cs = sqrt(float_t(8.0) * spade::consts::Rgas_uni * prim.T() * gas.mw_si[s] / spade::consts::pi);

				// Number density
				Num_den = spade::consts::Na_kmol * rhos[s] * gas.mw_si[s];

				// Relaxation time
				tau[s] = float_t(1.0) / (Num_den * cs * sigma);
			}
		}
		
		return tau;
	}

	// Computation for molar-averaged relaxation time
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_molarRelaxTime(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		using float_t = ptype;
		
		// Initialize relaxation time
		spade::ctrs::array<ptype, prim.nspecies()> tau = float_t(0.0);

		// Get molar concentrations
		spade::ctrs::array<ptype, prim.nspecies()> Xr = fluid_state::get_Xr(prim, gas);

		// Landau-teller interspecies relaxation time
		ptype sumXr, tau_sr, tau_sum;
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			// Check for molecules
			if (gas.isMol[s]>0)
			{
				sumXr   = float_t(0.0);
				tau_sum = float_t(0.0);
				// Sweep collision partners
				for (int r = 0; r<prim.nspecies(); ++r)
				{
					// Exclude electrons
					if (gas.mw_s[s]>1)
					{
						// Landau-teller relaxation time
						tau_sr = react.evaluate_tausr(s, r, prim.p(), prim.T(), gas);

						// Summation over r
						sumXr   += Xr[r];
						tau_sum += Xr[r] / tau_sr;
					}
				}

				// Compute relaxation time for species s
				tau[s] = sumXr / tau_sum;
			}
		}

		// Finished
		return tau;
	}
	
	// Computation for vibration energy exchange relaxation time
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_relaxationTime(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		// Initialize relaxation time components
		spade::ctrs::array<ptype, prim.nspecies()> tau_molar,tau_park;

		// Compute molar-averaged relaxation time
		tau_molar = compute_molarRelaxTime(prim, gas, react);

		// Compute Park relaxation time
		tau_park = compute_parkRelaxTime(prim, gas, react);

		// Total relaxation time
		return tau_molar + tau_park;
	}
	
	// Compute translational to vibrational source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	_sp_hybrid static ptype compute_St2v(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		// Vibrational energy
		spade::ctrs::array<ptype, prim.nspecies()> ev_st = get_evs(prim.T(), gas);
		spade::ctrs::array<ptype, prim.nspecies()> ev_s  = get_evs(prim.Tv(), gas);
			
		// Compute relaxation time
		spade::ctrs::array<ptype, prim.nspecies()> tau_s = compute_relaxationTime(prim, gas, react);

		// Compute species densities
		spade::ctrs::array<ptype, prim.nspecies()> rhos = fluid_state::get_rhos(prim, gas);

		// Compute source term component
		ptype St2v = 0.0;
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			if (gas.isMol[s]>0) St2v += rhos[s] * (ev_st[s] - ev_s[s]) / tau_s[s];
		}
		return St2v;
	}

	// Compute chemical to vibrational source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
    _sp_hybrid static ptype compute_Sc2v(const prim_chem_t<ptype, ns>& prim, const spade::ctrs::array<ptype, ns>& source, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		// Vibrational energy
		spade::ctrs::array<ptype, prim.nspecies()> ev_s  = get_evs(prim.Tv(), gas);
			
		// Compute source term component
		ptype Sc2v = 0.0;
		for (int s = 0; s<prim.nspecies(); ++s)
		{
			if (gas.isMol[s]>0) Sc2v += source[s] * ev_s[s];
		}
		return Sc2v;
	}

	// Compute heavy particle to electron source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
    _sp_hybrid static ptype compute_Sh2e(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		// Compute source term component
		return 0.0;
	}

	// Compute electron impact ionization source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib, const bool is_image>
	_sp_hybrid static ptype compute_Se2i(const prim_chem_t<ptype, ns>& prim, const spade::ctrs::array<ptype, ns>& source, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr, is_image>& react)
	{
		// Compute source term component
		return 0.0;
	}

	// Chemical source term structure. Primary driver
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib> struct chem_source_t
	{
		using info_type   = omni::info_list_t<omni::info::value>;
		using omni_type   = omni::prefab::cell_mono_t<info_type>;
		using output_type = fluid_state::flux_chem_t<ptype, ns>;
		using gas_t       = multicomponent_gas_t<ptype, ns, nvib>;
		using react_t     = reactionMechanism_t<ptype, ns, nr, true>;

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

			// Initialize source term
			spade::ctrs::array<ptype, q.nspecies()> source;
			ptype St2v, Sc2v, Sh2e, Se2i;

			// Compute chemical source term
			source = compute_chemSource(q, gas, reactions);
			
			// Compute vibrational source term components
			St2v = compute_St2v(q, gas, reactions);
			Sc2v = compute_Sc2v(q, source, gas, reactions);

			// Ionization-specific terms
			Sh2e = compute_Sh2e(q, gas, reactions);
			Se2i = compute_Se2i(q, source, gas, reactions);

			// Assign to out-going vector
			for (int s = 0; s<q.nspecies(); ++s) output.continuity(s) = source[s];
			output.energyVib() = St2v + Sc2v + Sh2e - Se2i;
			
			// Return vector for assignment onto the RHS
			return output;
		}
	};
	
}
