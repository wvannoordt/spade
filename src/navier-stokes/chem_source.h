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
	
	// Stores anything related to the employed reaction mechanism and/or vibrational relaxation model
	template<typename rtype, const std::size_t num_species, const std::size_t num_react> struct reactionMechanism_t
	{
		using float_t = rtype;
		// Constructor
		_sp_hybrid reactionMechanism_t(){}

		// Number of reactions in mechanism
		_sp_hybrid constexpr static std::size_t nspecies(){return num_species;}
		_sp_hybrid constexpr static std::size_t nreact(){return num_react;}

		// Reaction parameters
		constexpr static rtype Tmin=800; // Minimum reaction evaluation temperature
		constexpr static rtype eps_chem2=6400; // 80^2 --> Used in reaction temperature limiting
		constexpr static rtype sigma_coef=1E-20; // Collisional cross-section for Park relaxation time
		constexpr static rtype kb_max=1E26; // Maximum reaction rate
		spade::ctrs::array<int, num_species> nIntervals; // Over species
		spade::linear_algebra::dense_mat<rtype, num_species, 3> Tlower,Tupper; // Dynamic vectors storing temp intervals for gibbs curve fit
		spade::ctrs::array<int, num_react> reactType; // Over reactions
		spade::ctrs::array<rtype, num_react> Af,Ta,eta; // Over reactions
		spade::linear_algebra::dense_mat<rtype, num_species, 9> gibbsCoefT1,gibbsCoefT2,gibbsCoefT3; // Over species (9 coef per species)
		spade::linear_algebra::dense_mat<rtype, num_react, num_species> phi_diss=float_t(1.0); // Dissociation enhancement parameter
		spade::linear_algebra::dense_mat<rtype, num_react, num_species> vprod=float_t(0.0),vreact=float_t(0.0); // Stoichiometric coefficients
		
		// Vibrational relaxation parameters
		constexpr static rtype p0=101325.0;
		constexpr static rtype Asr=0.00116,Bsr=0.015; // For vibrational relaxation		
		
		// Forward controlling temperature
		_sp_hybrid rtype compute_Tf(const int& r, const rtype& T, const rtype& Tv) const
		{
			// Check reaction type
			if (reactType[r] == DISS_REACTION)
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
		_sp_hybrid rtype limiting(const rtype& T) const {return float_t(0.5) * ((T + Tmin) + sqrt(pow(T-Tmin,2) + eps_chem2));}
		
		// Compute arrhenius curve fit
		_sp_hybrid rtype compute_rates(const int& r, const rtype& T) const {return Af[r] * pow(T,eta[r]) * exp(-Ta[r] / T);}

		// Get necessary curve fitting coefficients for equilibrium constant
		_sp_hybrid void get_NASA9_coefficients(const int& s, const rtype& T, spade::ctrs::array<rtype, 9>& coefs) const
		{
			// Sweep temperature intervals
			for (int i = 0; i<nIntervals[s]; ++i)
			{
				// Check intervals
				if (T >= Tlower(s,i) && T < Tupper(s,i))
				{
					// Select interval
					if (i == 0)
					{
						// Interval 1
						for (int n = 0; n<coefs.size(); ++n) coefs[n] = gibbsCoefT1(s,n);
						return;
					}
					else if (i == 1)
					{
						// Interval 2
						for (int n = 0; n<coefs.size(); ++n) coefs[n] = gibbsCoefT2(s,n);
						return;
					}
					else if (i == 2)
					{
						// Interval 3
						for (int n = 0; n<coefs.size(); ++n) coefs[n] = gibbsCoefT3(s,n);
						return;
					}
				}
			}

			// Handle extremities
			if (T < Tlower(s,0))
			{
				// Interval 1
				for (int n = 0; n<coefs.size(); ++n) coefs[n] = gibbsCoefT1(s,n);
				return;
			}
			else if (T > Tupper(s,2))
			{
				// Interval 3
				for (int n = 0; n<coefs.size(); ++n) coefs[n] = gibbsCoefT3(s,n);
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
			rtype Acoef = Asr * sqrt(mu_sr) * pow(theta_min,float_t(4.0)/float_t(3.0));
			rtype Bcoef = Bsr * sqrt(sqrt(mu_sr));
			return (p0 / p) * exp(Acoef * (pow(T,-float_t(1.0)/float_t(3.0)) - Bcoef) - float_t(18.42));
		}
		
	};

	// Import gibbs energy data
	template<typename ptype, const std::size_t ns, const std::size_t nr>
	static void import_gibbsEnergy_data(const std::string& gibbs_fname, const std::vector<std::string>& speciesNames, reactionMechanism_t<ptype, ns, nr>& react)
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
				
				// Sweep species
				for (int s = 0; s<react.nspecies(); ++s)
				{
					// Do we need this species?
					if (species == speciesNames[s])
					{
						// Store lower/upper bound
						react.Tlower(s, nIntervals[s]) = Tlower;
						react.Tupper(s, nIntervals[s]) = Tupper;

						// Store coefficients
						if (nIntervals[s]==0)
						{
							for (int i = 0; i<gibbsCoef.size(); ++i) react.gibbsCoefT1(s,i) = gibbsCoef[i];
						}
						else if (nIntervals[s]==1)
						{
							for (int i = 0; i<gibbsCoef.size(); ++i) react.gibbsCoefT2(s,i) = gibbsCoef[i];
						}
						else if (nIntervals[s]==2)
						{
							for (int i = 0; i<gibbsCoef.size(); ++i) react.gibbsCoefT3(s,i) = gibbsCoef[i];
						}
						else
						{
							std::cerr << "Cannot hold more than 3 temperature intervals! Need to implement some 3D arrays maybe..." << std::endl;
						}
							
						// Count intervals
						nIntervals[s] += 1;
						
					}
				}
				// Check for end of file
				if (infile.eof()) break;
			}
			
			// Close file
			infile.close();
			
			// Copy intervals to reaction structure
			for (int s = 0; s<react.nspecies(); ++s) react.nIntervals[s] = nIntervals[s];
			
		}
		else
		{
			std::cerr << "Can not open provided gibbs energy data file!" << std::endl;
		}
		
		return;
	}

	// Import gibbs energy data
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	static void import_reaction_data(const std::string& react_fname, const std::string& gibbs_fname, const std::vector<std::string>& speciesNames, const multicomponent_gas_t<ptype, ns, nvib>& gas, reactionMechanism_t<ptype, ns, nr>& react)
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
					react.Af[numReact]  = Af;
					react.eta[numReact] = eta;
					react.Ta[numReact]  = Ta;

					// Set reaction type
					if (reactType == "diss")
					{
						react.reactType[numReact] = DISS_REACTION;
					}
					else if (reactType == "exch")
					{
						react.reactType[numReact] = EXCH_REACTION;
					}
					else if (reactType == "ei_ion")
					{
						react.reactType[numReact] = EI_ION_REACTION;
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
							if (participants[p] == speciesNames[s]) react.vreact(numReact, s) += float_t(1.0);
						}
					}
					
					// Set stoichiometric coefficients (products)
					for (int p = participants.size()/2; p<participants.size(); ++p)
					{
						for (int s = 0; s<react.nspecies(); ++s)
						{
							if (participants[p] == speciesNames[s]) react.vprod(numReact, s) += float_t(1.0);
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
								if (gas.isMol[s]>0) react.phi_diss(numReact,s) = phi;
							}
						}
						else if (species == "atoms")
						{
							// Enhance dissociation rates for atomic collisions
							for (int s = 0; s<react.nspecies(); ++s)
							{
								if (!gas.isMol[s]>0) react.phi_diss(numReact,s) = phi;
							}
						}
						else
						{
							// Enhance dissociation rates for specific collisions
							for (int s = 0; s<react.nspecies(); ++s)
							{
								if (species == speciesNames[s]) react.phi_diss(numReact,s) = phi;
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
				
				// Check for end of file
				if (infile.eof()) break;
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
		
		return;
	}

	// Compute forward reaction rates
	template<typename ptype, const std::size_t ns, const std::size_t nr>
	_sp_hybrid static void compute_forwardRates(const prim_chem_t<ptype, ns>& prim, const reactionMechanism_t<ptype, ns, nr>& react, spade::ctrs::array<ptype, nr>& kf, spade::ctrs::array<ptype, nr>& kfb)
	{
		// Loop reactions
		ptype Tf, Tb;
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
			kfb[r] = react.compute_rates(r, Tb);
		}

		return;
	}

	// Compute forward reaction rates
	template<typename ptype, const std::size_t ns, const std::size_t nr>
	_sp_hybrid static spade::ctrs::array<ptype, nr> compute_backwardRates(const prim_chem_t<ptype, ns>& prim, const reactionMechanism_t<ptype, ns, nr>& react, const spade::ctrs::array<ptype, nr>& kfb)
	{
		using float_t = ptype;

		// Initialize
		spade::ctrs::array<ptype, react.nreact()> kb = 0.0;
		ptype Tb,hi,si,gibbs,Kc,vr;
		spade::ctrs::array<ptype, 9> coefs;
		ptype T,T2,T3,T4,logT;
		
		// Sweep reactions
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Initialize gibbs energy
			gibbs = float_t(0.0);

			// Initialize stoichiometric summation
			vr = float_t(0.0);

			// Backward controlling temperature
			Tb = react.compute_Tb(r, prim.T(), prim.Tv());
			Tb = react.limiting(Tb);

			// Pre-compute some values
			T    = Tb;
			T2   = Tb*Tb;
			T3   = T2*Tb;
			T4   = T3*Tb;
			logT = log(Tb);
			
			// Sweep species
			for (int s = 0; s<prim.nspecies(); ++s)
			{
				// Get curve fitting coefficients
				react.get_NASA9_coefficients(s, Tb, coefs);
				
				// Enthalpy norm
				hi = react.evaluate_enthalpyNorm(s, Tb, coefs);
				
				// Entropy norm
				si = react.evaluate_entropyNorm(s, Tb, coefs);
				
				// Gibbs energy summation
				gibbs += (react.vprod(r,s) - react.vreact(r,s)) * (hi - si);
				
				// Stoichiometric summation
				vr += react.vprod(r,s) - react.vreact(r,s);
			}

			// Compute equilibrium constant
			Kc = exp(-gibbs) * pow(float_t(0.1) / (spade::consts::Rgas_uni * Tb * 1E-3),vr);

			// Compute backward reaction rate
			const auto kb_max = react.kb_max;
			kb[r] = utils::min(kfb[r] / Kc, kb_max);
		}

		return kb;
	}

	// Compute reaction product
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_reactionProduct(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react, const spade::ctrs::array<ptype, nr>& kf, const spade::ctrs::array<ptype, nr>& kb)
	{
		using float_t = ptype;
		// Initialize some variables
		spade::ctrs::array<ptype, react.nreact()> con   = float_t(0.0);
		spade::ctrs::array<ptype, react.nreact()> Rf_kf = float_t(1000.0);
		spade::ctrs::array<ptype, react.nreact()> Rb_kb = float_t(1000.0);

		// Initialize source term
		spade::ctrs::array<ptype, prim.nspecies()> source = float_t(0.0);

		// Compute species density
		spade::ctrs::array<ptype, prim.nspecies()> rhos = fluid_state::get_rhos(prim, gas);

		// Sweep reactions
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Check for dissociation reactions
			if (react.reactType[r] == DISS_REACTION)
			{
				// Sweep species
				for (int s = 0; s<prim.nspecies(); ++s)
				{
					// Block out electrons
					if (gas.mw_s[s]>1) con[r] += rhos[s] * react.phi_diss(r,s) * gas.mw_si[s] * 1E-3;
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
				if (react.vreact(r,s) != float_t(0.0)) Rf_kf[r] *= pow(1E-3 * rhos[s] * gas.mw_si[s], react.vreact(r,s));

				// Backward reaction rate
				if (react.vprod(r,s) != float_t(0.0)) Rb_kb[r] *= pow(1E-3 * rhos[s] * gas.mw_si[s], react.vprod(r,s));
			}
		}

		// Final computation for source term
		for (int r = 0; r<react.nreact(); ++r)
		{
			// Sweep species in reaction
			for (int s = 0; s<prim.nspecies(); ++s)
			{
				source[s] += (react.vprod(r,s) - react.vreact(r,s)) * (kf[r] * Rf_kf[r] - kb[r] * Rb_kb[r]) * gas.mw_s[s] * con[r];
			}
		}

		return source;
	}

	// Compute chemical source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_chemSource(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
	{
		// Allocate vectors for reactions
		spade::ctrs::array<ptype, react.nreact()> kf,kfb,kb;
		spade::ctrs::array<ptype, prim.nspecies()> source;
		
		// Compute forward reaction rates
		compute_forwardRates(prim, react, kf, kfb);
		
		// Compute backward reaction rates
		kb = compute_backwardRates(prim, react, kfb);

		// Compute reaction product
		source = compute_reactionProduct(prim, gas, react, kf, kb);
		
		return source;
	}

	// Computation for Park vibrational relaxation time for high temperature corrections
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_parkRelaxTime(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
	{
		using float_t = ptype;
		
		// Initialize relaxation time
		spade::ctrs::array<ptype, prim.nspecies()> tau = float_t(0.0);
		
		// Compute collision cross-section
		ptype sigma = react.sigma_coef * float_t(2.5E9) / (prim.T() * prim.T());

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
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_molarRelaxTime(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
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
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static spade::ctrs::array<ptype, ns> compute_relaxationTime(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
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
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static ptype compute_St2v(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
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
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
    _sp_hybrid static ptype compute_Sc2v(const prim_chem_t<ptype, ns>& prim, const spade::ctrs::array<ptype, ns>& source, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
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
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
    _sp_hybrid static ptype compute_Sh2e(const prim_chem_t<ptype, ns>& prim, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
	{
		// Compute source term component
		return 0.0;
	}

	// Compute electron impact ionization source term
	template<typename ptype, const std::size_t ns, const std::size_t nr, const std::size_t nvib>
	_sp_hybrid static ptype compute_Se2i(const prim_chem_t<ptype, ns>& prim, const spade::ctrs::array<ptype, ns>& source, const multicomponent_gas_t<ptype, ns, nvib>& gas, const reactionMechanism_t<ptype, ns, nr>& react)
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
		using react_t     = reactionMechanism_t<ptype, ns, nr>;

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
