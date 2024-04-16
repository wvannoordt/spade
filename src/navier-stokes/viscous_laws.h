#pragma once

#include "core/c20.h"
#include "omni/omni.h"
#include "navier-stokes/fluid_state.h"


namespace spade::viscous_laws
{
    template <class T> concept state_independent_viscosity = std::floating_point<typename T::value_type> && requires(T t)
    {
        typename T::value_type;
        t.get_visc();          // mu
        t.get_beta();          // -2/3 * mu
        t.get_diffuse();       // cp * mu/Pr
    };
    
    template <class T> concept state_dependent_viscosity = std::floating_point<typename T::value_type>
    && requires(T t, const fluid_state::prim_t<typename T::value_type>& q)
    {
        typename T::value_type;
        t.get_visc(q);
        t.get_beta(q);
        t.get_diffuse(q);
    };
    
    template <class T> concept viscous_law = state_dependent_viscosity<T> || state_independent_viscosity<T>;

    template <typename data_t>
    struct visc_result_t
    {
        data_t mu;    // viscosity
        data_t beta;  // second viscosity
        data_t alpha; // heat conductivity
    };
    
    template <typename derived_t, typename ilist_t = omni::info_list_t<>>
    struct visc_law_interface_t
    {
        _sp_hybrid derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        _sp_hybrid const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type   = ilist_t;
        // using result_type = visc_result_t<typename derived_t::value_type>;

        _sp_hybrid auto get_visc        (const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_visc        (args...);}, input);
        }

        _sp_hybrid auto get_beta        (const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_beta        (args...);}, input);
        }

        _sp_hybrid auto get_diffuse(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_diffuse(args...);}, input);
        }
    };
    
    template <typename dtype> struct constant_viscosity_t
    : public visc_law_interface_t<constant_viscosity_t<dtype>>
    {
        using base_t = visc_law_interface_t<constant_viscosity_t<dtype>>;
        using value_type = dtype;
        using result_type = visc_result_t<value_type>;
        using base_t::get_visc;
        using base_t::get_beta;
        using base_t::get_diffuse;
        using base_t::info_type;

        dtype visc;
        dtype prandtl_inv;
        dtype beta;
        
        constant_viscosity_t(const dtype& visc_in, const dtype& prandtl_in)
        : visc{visc_in}, beta{dtype(-2.0*visc_in/3.0)}, prandtl_inv{dtype(1.0)/prandtl_in}
        { }
        
        _sp_hybrid result_type get_all(const auto&) const
        {
            return result_type{this->get_visc(), this->get_beta(), this->get_diffuse()};
        }
        
        _sp_hybrid dtype get_visc() const
        {
            return visc;
        }
        
        _sp_hybrid dtype get_beta() const
        {
            return beta;
        }
        
        _sp_hybrid dtype get_diffuse() const
        {
            return visc*prandtl_inv;
        }
    };
    
    template <typename dtype> struct power_law_t
    : public visc_law_interface_t<power_law_t<dtype>, omni::info_list_t<omni::info::value>>
    {
            typedef dtype value_type;
            dtype mu_ref;
            dtype T_ref;
            dtype power;
            dtype prandtl;

            using base_t = visc_law_interface_t<power_law_t<dtype>, omni::info_list_t<omni::info::value>>;
            using base_t::get_visc;
            using base_t::get_beta;
            using base_t::get_diffuse;
            using base_t::info_type;
            
            power_law_t(const dtype& mu_ref_in, const dtype& T_ref_in, const dtype& power_in, const dtype& prandtl_in)
            : mu_ref{mu_ref_in}, T_ref{T_ref_in}, power{power_in}, prandtl{prandtl_in}
            { }
            
            template <fluid_state::is_state_type state_t> _sp_hybrid dtype get_visc(const state_t& q) const
            {
                return mu_ref*std::pow(q.T()/T_ref, power);
            }
            
            template <fluid_state::is_state_type state_t> _sp_hybrid dtype get_beta(const state_t& q) const
            {
                return -0.66666666667*this->get_visc(q);
            }
            
            template <fluid_state::is_state_type state_t> _sp_hybrid dtype get_diffuse(const state_t& q) const
            {
                return this->get_visc(q)/prandtl;
            }
    };
    
    template <typename state_t, typename visc_func_t, typename beta_func_t, typename cond_func_t> struct udf_t
    : public visc_law_interface_t<udf_t<state_t, visc_func_t, beta_func_t, cond_func_t>, omni::info_list_t<omni::info::value>>
    {
            typedef typename state_t::value_type value_type;

            using base_t = visc_law_interface_t<udf_t<state_t, visc_func_t, beta_func_t, cond_func_t>, omni::info_list_t<omni::info::value>>;
            using base_t::get_visc;
            using base_t::get_beta;
            using base_t::get_diffuse;
            using base_t::info_type;
            
            const visc_func_t& vfunc;
            const beta_func_t& bfunc;
            const cond_func_t& cfunc;
            
            udf_t(const state_t& state, const visc_func_t& v_in, const beta_func_t& b_in, const cond_func_t& c_in)
            : vfunc{v_in},
            bfunc{b_in},
            cfunc{c_in}
            {}
            
            _sp_hybrid value_type get_visc(const state_t& q) const
            {
                return vfunc(q);
            }
            
            _sp_hybrid value_type get_beta(const state_t& q) const
            {
                return bfunc(q);
            }
            
            _sp_hybrid value_type get_diffuse(const state_t& q) const
            {
                return cfunc(q);
            }
    };

    template <typename laminar_t, typename turb_t> struct sgs_visc_t
    {
        using info_type = omni::info_union<typename laminar_t::info_type, typename turb_t::info_type>;
        using value_type = decltype(typename laminar_t::value_type() + typename turb_t::value_type());
        using result_type = visc_result_t<value_type>;

        const laminar_t lam;
        const turb_t turb;
        sgs_visc_t(const laminar_t& lam_in, const turb_t& turb_in) : lam{lam_in}, turb{turb_in} {}
        
        _sp_hybrid result_type get_all(const auto& data) const
        {
            const auto lam_all  = lam.get_all(data);
            
            const auto mu_t     = turb.get_mu_t(data);
            const auto pr_t     = turb.get_prt(data);
            
            return result_type{
                lam_all.mu    + mu_t,
                lam_all.beta  - value_type(0.66666666667)*mu_t,
                lam_all.alpha + mu_t/pr_t};
        }

        _sp_hybrid value_type get_visc(const auto& data) const
        {
            const auto lv = lam.get_visc(data);
            const auto tv = turb.get_mu_t(data);
            return lv + tv;
        }
        
        _sp_hybrid value_type get_beta(const auto& data) const
        {
            return -0.66666666667*this->get_visc(data);
        }
        
        _sp_hybrid value_type get_diffuse(const auto& data) const
        {
            const auto ld = lam.get_diffuse(data);
            const auto td = turb.get_mu_t(data)/turb.get_prt(data);
            return ld + td;
        }
    };

    // Gupta Viscous Model (Multicomponent Gas) <JRB | Implemented: 4-14-24 | Validated: TODO>
    template <typename dtype, const std::size_t ns, const std::size_t max_vib> struct gupta_visc_t
    : public visc_law_interface_t<gupta_visc_t<dtype, ns, max_vib>, omni::info_list_t<omni::info::value>>
    {
        typedef dtype value_type;
        typedef fluid_state::multicomponent_gas_t<dtype, ns, max_vib> gas_type;

        using base_t = visc_law_interface_t<gupta_visc_t<dtype, ns, max_vib>, omni::info_list_t<omni::info::value>>;
        using base_t::get_visc;
        using base_t::get_beta;
        using base_t::get_diffuse;
        using base_t::info_type;
        
        _sp_hybrid constexpr static std::size_t num_species(){return ns;}
        _sp_hybrid constexpr static std::size_t num_collisions(){return std::size_t((ns*(ns+1))/2);}
        _sp_hybrid constexpr static std::size_t map_indices(const std::size_t s1, const std::size_t s2)
        {
            return std::size_t(ns*s1 + s2 - (s1*(s1+1))/2);
        }

        // the coulomb collision fit parameters (can't think of a cleaner way to store these)
        constexpr static dtype c1_att = 0.0313, C1_att = -0.476, D1_att = 0.784;
        constexpr static dtype c1_rep = 0.0106, C1_rep = 0.138,  D1_rep = 0.765;
        constexpr static dtype c2_att = 0.0377, C2_att = -0.146, D2_att = 1.262;
        constexpr static dtype c2_rep = 0.0274, C2_rep = 0.157,  D2_rep = 1.235;

        // the fit parameter arrays, with n(n+1)/2 elements, corresponding to each unique combination of species
        spade::ctrs::array<dtype, std::size_t((ns*(ns+1))/2)> A0, A1, A2, A3;

        const fluid_state::multicomponent_gas_t<dtype, ns, max_vib>& gas;

        // store the species molecular weight and charge
        gupta_visc_t(const gas_type& gas_in) : gas{gas_in} { }
        
        _sp_hybrid dtype get_visc(const fluid_state::prim_chem_t<dtype, ns>& q) const
        {
            // the output variable
            dtype mu_out = 0;

            // get the molar concentrations and density
            const spade::ctrs::array<dtype, ns> Xr = fluid_state::get_Xr(q, gas);
            const spade::ctrs::array<dtype, ns> rhos = fluid_state::get_rhos(q, gas);

            // get the electron number density in (#/cm^3)
            dtype n_e = 0;
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                if (gas.mw_s[s1] < 1)
                {
                    n_e = std::pow(10, 6)*consts::Na_kmol*rhos[s1]/gas.mw_s[s1];
                    break;
                }
            }

            for (int s1 = 0; s1 < num_species(); s1++)
            {
                // species molecular mass in (kg)
                const dtype mass_s  = gas.mw_s[s1]/consts::Na_kmol;         

                // the denominator term
                dtype denominator = 0;
                
                // compute the collision terms at different controlling temperatures
                for (int s2 = 0; s2 < num_species(); s2++)
                {
                    dtype delta2 = 0;
                    dtype omega2 = 0;

                    if (gas.mw_s[s1] < 1)
                    {
                        // an electron
                        if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                        {
                            // neutral-electron (vibrational/electronic controlling temperature)
                            omega2 = A3[map_indices(s1, s2)]*std::pow(q.Tv(), A0[map_indices(s1, s2)]*std::log(q.Tv())*std::log(q.Tv()) + A1[map_indices(s1, s2)]*std::log(q.Tv()) + A2[map_indices(s1, s2)]);
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.Tv()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                        else
                        {
                            // ion-electron or electron-electron (vibrational/electronic controlling temperature)
                            const dtype lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                            const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));

                            dtype c2, C2, D2;
                            if (spade::utils::sign(gas.charge_s[s1]) != spade::utils::sign(gas.charge_s[s2])) 
                            { c2 = c2_att; C2 = C2_att; D2 = D2_att; }
                            else 
                            { c2 = c2_rep; C2 = C2_rep; D2 = D2_rep; }

                            omega2 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.Tv()*q.Tv()))*std::log(D2*T_star*(dtype(1)-C2*std::exp(-1*c2*T_star))+dtype(1));
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.Tv()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                    }
                    else 
                    {
                        // not electron (translational controlling temperature)
                        if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                        {
                            // neutral
                            omega2 = A3[map_indices(s1, s2)]*std::pow(q.T(), A0[map_indices(s1, s2)]*std::log(q.T())*std::log(q.T()) + A1[map_indices(s1, s2)]*std::log(q.T()) + A2[map_indices(s1, s2)]);
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.T()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                        else if (gas.mw_s[s2] < 1)
                        {
                            // ion-electron (vibrational/electronic controlling temperature)
                            const dtype lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                            const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));

                            dtype c2, C2, D2;
                            if (spade::utils::sign(gas.charge_s[s1]) != spade::utils::sign(gas.charge_s[s2])) 
                            { c2 = c2_att; C2 = C2_att; D2 = D2_att; }
                            else 
                            { c2 = c2_rep; C2 = C2_rep; D2 = D2_rep; }

                            omega2 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.Tv()*q.Tv()))*std::log(D2*T_star*(dtype(1)-C2*std::exp(-1*c2*T_star))+dtype(1));
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.Tv()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                        else
                        {
                            // ion-ion (translational controlling temperature)
                            const dtype lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                            const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));

                            dtype c2, C2, D2;
                            if (spade::utils::sign(gas.charge_s[s1]) != spade::utils::sign(gas.charge_s[s2])) 
                            { c2 = c2_att; C2 = C2_att; D2 = D2_att; }
                            else 
                            { c2 = c2_rep; C2 = C2_rep; D2 = D2_rep; }
                            
                            omega2 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.T()*q.T()))*std::log(D2*T_star*(dtype(1)-C2*std::exp(-1*c2*T_star))+dtype(1));
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.T()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                    }

                    // add the result to the denominator term
                    denominator += Xr[s2]*delta2;
                }

                // now add the result to the viscosity
                mu_out += mass_s*Xr[s1] / denominator;
            }

            return mu_out;
        }
        
        _sp_hybrid dtype get_beta(const fluid_state::prim_chem_t<dtype, ns>& q) const
        {
            return -0.66666666667*this->get_visc(q);
        }
        
        _sp_hybrid spade::ctrs::array<dtype, ns> get_diffuse(const fluid_state::prim_chem_t<dtype, ns>& q) const
        {
            spade::ctrs::array<dtype, ns> diffuse_out;

            // get the molar concentrations and density
            const spade::ctrs::array<dtype, ns> Xr = fluid_state::get_Xr(q, gas);
            const spade::ctrs::array<dtype, ns> rhos = fluid_state::get_rhos(q, gas);

            // compute the sum of molar concentrations
            dtype Xr_sum = 0;
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                Xr_sum += Xr[s1];
            }

            // get the electron number density in (#/cm^3)
            dtype n_e = 0;
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                if (gas.mw_s[s1] < 1)
                {
                    n_e = std::pow(10, 6)*consts::Na_kmol*rhos[s1]/gas.mw_s[s1];
                    break;
                }
            }

            for (int s1 = 0; s1 < num_species(); s1++)
            {
                // the denominator term
                dtype denominator = 0;
                
                // compute the collision terms at different controlling temperatures
                for (int s2 = 0; s2 < num_species(); s2++)
                {
                    if (s2 != s1)
                    {
                        // binary diffusion coefficient
                        dtype binary_diff = 0;
                        dtype delta1      = 0;
                        dtype omega1      = 0;

                        // if the collision involves an electron (use vib/elec temperature)
                        if (gas.mw_s[s1] < 1 || gas.mw_s[s2] < 1)
                        {
                            if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                            {
                                // neutral+electron (use fit params)
                                omega1 = A3[map_indices(s1, s2)]*std::pow(q.Tv(), A0[map_indices(s1, s2)]*std::log(q.Tv())*std::log(q.Tv()) + A1[map_indices(s1, s2)]*std::log(q.Tv()) + A2[map_indices(s1, s2)]);
                            }
                            else
                            {
                                // electron+ion (use colomb params)
                                const dtype lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                                const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));
                                
                                dtype c1, C1, D1;
                                if (spade::utils::sign(gas.charge_s[s1]) != spade::utils::sign(gas.charge_s[s2])) 
                                { c1 = c1_att; C1 = C1_att; D1 = D1_att; }
                                else 
                                { c1 = c1_rep; C1 = C1_rep; D1 = D1_rep; }
                            
                                omega1 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.Tv()*q.Tv()))*std::log(D1*T_star*(dtype(1)-C1*std::exp(-1*c1*T_star))+dtype(1));
                            }
                            delta1 = dtype(8)/3*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/ (consts::pi*consts::Rgas_uni*q.Tv()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10, -20)*omega1;
                            binary_diff = consts::kSI*q.Tv()/(q.p()*delta1);
                        }
                        // otherwise (use trans/rot temperature)
                        else 
                        {
                            if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                            {
                                // neutral+neutral/ion (use fit params)
                                omega1 = A3[map_indices(s1, s2)]*std::pow(q.T(), A0[map_indices(s1, s2)]*std::log(q.T())*std::log(q.T()) + A1[map_indices(s1, s2)]*std::log(q.T()) + A2[map_indices(s1, s2)]);
                            }
                            else
                            {
                                // ion+ion (use colomb params)
                                const dtype lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                                const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));

                                dtype c1, C1, D1;
                                if (spade::utils::sign(gas.charge_s[s1]) != spade::utils::sign(gas.charge_s[s2])) 
                                { c1 = c1_att; C1 = C1_att; D1 = D1_att; }
                                else 
                                { c1 = c1_rep; C1 = C1_rep; D1 = D1_rep; }
                            
                                omega1 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.T()*q.T()))*std::log(D1*T_star*(dtype(1)-C1*std::exp(-1*c1*T_star))+dtype(1));
                            }
                            delta1 = dtype(8)/3*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/ (consts::pi*consts::Rgas_uni*q.T()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10, -20)*omega1;
                            binary_diff = consts::kSI*q.T()/(q.p()*delta1);
                        }

                        // add the term to the denominator
                        denominator += Xr[s2]/binary_diff;
                    }
                }

                // now compute and store the species diffusion coefficient
                diffuse_out[s1] = Xr_sum*Xr_sum*gas.mw_s[s1]*(1-gas.mw_s[s1]*Xr[s1])/denominator;
            }

            return diffuse_out;
        }
    };

    // Initialization Function for the Collision Integral Fit Parameters <JRB | Implemented: 4-14-24 | Validated: TODO>
	template<typename dtype, const std::size_t ns, const std::size_t max_vib>
	static void import_gupta_collision_integral_data(const std::string& fname, const std::vector<std::string>& speciesNames, fluid_state::multicomponent_gas_t<dtype, ns, max_vib>& gas, gupta_visc_t<dtype, ns>& gupta_visc)
	{
		// Open the input file
        std::ifstream infile;
		std::string full_fname = std::getenv("SPADE") + std::string("/src/navier-stokes/speciesCollision/") + fname;
        infile.open(full_fname);
		
		if (infile)
		{
			// variables to temporarily store data
			std::string species1, species2;
            dtype A0, A1, A2, A3;
            spade::ctrs::array<std::size_t, ns*ns> unique_indices, unique_indices_check;

			// sweep entire file
			while (true)
			{
				// read line in file
				infile >> species1 >> species2 >> A0 >> A1 >> A2 >> A3;
				
				// sweep through each combination of species
				for (int s1 = 0; s1 < speciesNames.size(); s1++)
				{
                    for (int s2 = 0; s2 < speciesNames.size(); s2++)
				    {
                        // Check to see if either permutation is found, and assign the data to the corresponding index
                        if (species1 == speciesNames[s1] && species2 == speciesNames[s2])
                        {
                            // also check to make sure the collision is neutral
                            if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                            {
                                gupta_visc.A0[gupta_visc.map_indices(s1, s2)] = A0;
                                gupta_visc.A1[gupta_visc.map_indices(s1, s2)] = A1;
                                gupta_visc.A2[gupta_visc.map_indices(s1, s2)] = A2;
                                gupta_visc.A3[gupta_visc.map_indices(s1, s2)] = A3;

                                unique_indices[s1*ns+s2] = gupta_visc.map_indices(s1, s2);
                                unique_indices[s2*ns+s1] = gupta_visc.map_indices(s1, s2);
                            }
                        }
                    }
				}

                // check for end of file
				if (infile.eof()) break;
			}
			
			// close file
			infile.close();

            // check to make sure all neutral collisions are loaded in
            for (int s1 = 0; s1 < speciesNames.size(); s1++)
            {
                for (int s2 = 0; s2 < speciesNames.size(); s2++)
                {
                    if(gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                    {
                        unique_indices_check[s1*ns+s2] = gupta_visc.map_indices(s1, s2);
                        unique_indices_check[s2*ns+s1] = gupta_visc.map_indices(s1, s2);
                    }
                }
            }
            if (unique_indices != unique_indices_check)
            {
                std::cerr << "Missing neutral collision data in collision file!" << std::endl;
            }
		}
		else
		{
			std::cerr << "Can not open provided species collision file!" << std::endl;
		}

		return;
	}

}