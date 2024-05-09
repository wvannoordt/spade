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

    // WIP CODE
    //
    //

    // actually define the different fit types up here (as structs), and have them contain the relavent info
    
    /*
    template <typename dtype, const std::size_t num_fit_params> struct gupta_collision_fit_type_t
    {
        _sp_hybrid constexpr static std::size_t num_params() { return num_fit_params; }
        _sp_hybrid virtual constexpr static dtype fit_func() = 0;



    }*/

    // TEMPORARY NAMES <JRB>
    enum gupta_collision_fit_type
    {
        coulomb    = 0,
        pow_log2   = 1,
        frac_T_T2  = 2,
        frac_T2_T3 = 3,
        exp_log3   = 4,
    }; 

    // overloading the stream operator for this enumerated type
    std::istream& operator >> (std::istream& in, gupta_collision_fit_type& fit_type) 
    {
        std::string value;
        if(in >> value)
        {
            if      (value == "0") fit_type = coulomb;
            else if (value == "1") fit_type = pow_log2;
            else if (value == "2") fit_type = frac_T_T2;
            else if (value == "3") fit_type = frac_T2_T3;
            else if (value == "4") fit_type = exp_log3;
        }
        return in;
    };

    // Gupta Collision Integral Fit Functions <JRB>
    template <typename dtype, typename gas_t> struct gupta_collision_fit_t
    {
        typedef dtype float_t;

        // the maximum number of collision integral fit parameters
        constexpr static std::size_t max_num_params = 6;
        
        _sp_hybrid constexpr static std::size_t num_species()    { return gas_t::nspecies(); }
        //_sp_hybrid constexpr static std::size_t num_params()     { return num_fit_params; }
        _sp_hybrid constexpr static std::size_t num_collisions() { return std::size_t((num_species()*(num_species()+1))/2); }

        // this function maps down a triangular matrix into a linear array
        _sp_hybrid constexpr static std::size_t map_indices(const std::size_t s1, const std::size_t s2)
        {
            return std::size_t(num_species()*std::min(s1, s2) + std::max(s2, s1) - (std::min(s1, s2)*(std::min(s1, s2)+1))/2);
        }

        // this function returns the number of parameters used by different fit types
        _sp_hybrid constexpr static std::size_t num_params(const gupta_collision_fit_type fit_type)
        {
            switch (fit_type) 
            {
                case coulomb:    return 3;
                case pow_log2:   return 4;
                case frac_T_T2:  return 4;
                case frac_T2_T3: return 6;
                case exp_log3:   return 4;
            }
            return 0;
        }

        // the fit parameter matrices, with n(n+1)/2 elements, corresponding to each unique combination of species
        spade::linear_algebra::dense_mat<float_t, max_num_params, std::size_t((num_species()*(num_species()+1))/2)> omega11_params, omega22_params;
        gupta_collision_fit_type fit_types[std::size_t((num_species()*(num_species()+1))/2)];

        // the coulomb shieliding collision fit parameters (not sure if there's a better way of storing these)
        // float_t coulomb_fit_params_11[2][3], coulomb_fit_params_11[2][3];

        // constexpr static float_t coulomb_fit_params_11[2][3] = {{0.0313, -0.476, 0.784}, {0.0106, 0.138, 0.765}}; // {{att}, {rep}}
        // constexpr static float_t coulomb_fit_params_22[2][3] = {{0.0377, -0.146, 1.262}, {0.0274, 0.157, 1.235}}; // {{att}, {rep}}

        // the constructor
        gupta_collision_fit_t() 
        {
            // the fit type and number of fit parameters need to be consistent
            //if constexpr(fit_type == pow_log2 || fit_type == exp_log3 || fit_type == frac_T_T2) 
            //    static_assert(num_fit_params == 4, "Selected collision integral fit type requires four fit parameters");
            //else if constexpr(fit_type == frac_T2_T3) 
            //    static_assert(num_fit_params == 6, "Selected collision integral fit type requires six fit parameters"); 
        }

        // this loads in the fit data from the specified data files
        _sp_hybrid void load_collision_integral_data(const std::string& neutral_fname, const std::string& ion_fname, const std::vector<std::string>& speciesNames, gas_t& gas)
        {
            // Open the neutrals input file
            std::ifstream neutrals_file;
            std::string full_neutral_fname = std::getenv("SPADE") + std::string("/src/navier-stokes/speciesTransport/") + neutral_fname; // we need a better way of accessing data files
            neutrals_file.open(full_neutral_fname);
            constexpr std::size_t skip_neutral_lines = 1;

            // Open the ions input file
            std::ifstream ions_file;
            std::string full_ion_fname = std::getenv("SPADE") + std::string("/src/navier-stokes/speciesTransport/") + ion_fname; // we need a better way of accessing data files
            ions_file.open(full_ion_fname);
            constexpr std::size_t skip_ion_lines = 2;
            
            // read in the collision integral fit data from the neutrals input file
            if (neutrals_file)
            {
                // skip the specified number of lines
                for (int i = 0; i < skip_neutral_lines; i++)
                    neutrals_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // variables to temporarily store data
                std::string species1, species2;
                gupta_collision_fit_type fit_type;
                std::string integral_type;
                float_t collision_params[max_num_params];

                // check arrays to make sure each species collision is accounted for
                spade::ctrs::array<std::size_t, num_species()*num_species()> unique_indices, unique_indices_check;

                // sweep entire file
                while (true)
                {
                    // collision file format
                    neutrals_file >> species1 >> species2 >> fit_type >> integral_type;
                    for (int i = 0; i < num_params(fit_type); i++)
                        neutrals_file >> collision_params[i];
                    
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
                                    fit_types[map_indices(s1, s2)] = fit_type;

                                    for (int i = 0; i < num_params(fit_type); i++)
                                    {
                                        if (integral_type == "11")
                                            omega11_params(i, map_indices(s1, s2)) = collision_params[i];
                                        else if (integral_type == "22")
                                            omega22_params(i, map_indices(s1, s2)) = collision_params[i];
                                    }

                                    unique_indices[s1*num_species()+s2] = map_indices(s1, s2);
                                    unique_indices[s2*num_species()+s1] = map_indices(s1, s2);
                                }
                            }
                        }
                    }

                    // check for end of file
                    if (neutrals_file.eof()) break;
                }
                
                // close neutrals file
                neutrals_file.close();

                // check to make sure all neutral collisions are loaded in
                for (int s1 = 0; s1 < speciesNames.size(); s1++)
                {
                    for (int s2 = 0; s2 < speciesNames.size(); s2++)
                    {
                        if(gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                        {
                            unique_indices_check[s1*num_species()+s2] = map_indices(s1, s2);
                            unique_indices_check[s2*num_species()+s1] = map_indices(s1, s2);
                        }
                    }
                }
                /*if (unique_indices != unique_indices_check)
                    std::cout << "Missing neutral collision data in collision file!" << std::endl;*/
            }
            /*else std::cout << "Can not open provided neutral species collision file!" << std::endl;*/ // NEED TO THROW EXCEPTIONS HERE?

            // read in the data from the ions file
            if (ions_file)
            {
                // skip the specified number of lines
                for (int i = 0; i < skip_ion_lines; i++)
                    ions_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // temporary variables to parse the data
                std::string integral_type;
                float_t coulomb_fit_params[2][3];

                // sweep entire file
                while (true)
                {
                    // now read in the shielding fit parameters
                    ions_file >> integral_type;
                    for (int i = 0; i < num_params(coulomb); i++)
                        for (int j = 0; j < 2; j++)
                            ions_file >> coulomb_fit_params[j][i];
                    
                    // sweep through each combination of species
                    for (int s1 = 0; s1 < speciesNames.size(); s1++)
                    {
                        for (int s2 = 0; s2 < speciesNames.size(); s2++)
                        {
                            // check to see if the collision is not neutral
                            if (gas.charge_s[s1] != 0 && gas.charge_s[s2] != 0)
                            {
                                // add the fit type to the corresponding array
                                fit_types[map_indices(s1, s2)] = coulomb;

                                // check if the collision is attractive or repulsive
                                int rep = 0;
                                if (gas.charge_s[s1] == gas.charge_s[s2]) 
                                    rep = 1;

                                for (int i = 0; i < num_params(coulomb); i++)
                                {
                                    if (integral_type == "11")
                                        omega11_params(i, map_indices(s1, s2)) = coulomb_fit_params[rep][i];
                                    else if (integral_type == "22")
                                        omega22_params(i, map_indices(s1, s2)) = coulomb_fit_params[rep][i];
                                }
                            }
                        }
                    }

                    // check for end of file
                    if (ions_file.eof()) break;
                }
                
                // close neutrals file
                ions_file.close();
            }
            //else std::cout << "Can not open provided ion collision file!" << '\n';



        }

        // collision integral for neutrals
        _sp_hybrid float_t omega_fit(const std::size_t& s1, const std::size_t& s2, const float_t& T, const auto& omega_params) const
        {
            const gupta_collision_fit_type fit_type = fit_types[map_indices(s1, s2)];

            float_t omega_nn;
            if (fit_type == pow_log2)
                omega_nn = omega_params(3, map_indices(s1, s2))*std::pow(T, omega_params(0, map_indices(s1, s2))*std::log(T)*std::log(T) + omega_params(1, map_indices(s1, s2))*std::log(T) + omega_params(2, map_indices(s1, s2)));
            else if (fit_type == exp_log3)
                omega_nn = std::exp(omega_params(0, map_indices(s1, s2))*std::log(T)*std::log(T)*std::log(T) + omega_params(2, map_indices(s1, s2))*std::log(T)*std::log(T) + omega_params(3, map_indices(s1, s2))*std::log(T) + omega_params(0, map_indices(s1, s2))*std::log(T));
            else if (fit_type == frac_T_T2)
                omega_nn = (omega_params(0, map_indices(s1, s2))*T + omega_params(1, map_indices(s1, s2))) / (T*T + omega_params(2, map_indices(s1, s2))*T + omega_params(3, map_indices(s1, s2)));
            else if (fit_type == frac_T2_T3)
                omega_nn = (omega_params(0, map_indices(s1, s2))*T*T + omega_params(1, map_indices(s1, s2))*T + omega_params(2, map_indices(s1, s2))) / (T*T*T + omega_params(3, map_indices(s1, s2))*T*T + omega_params(4, map_indices(s1, s2))*T + omega_params(5, map_indices(s1, s2)));
            else omega_nn = 0.0; // FIT NOT IMPLEMENTED

            return omega_nn;
        }

        // collision integral using coulomb shielding fit for ion-ion and ion-electron
        _sp_hybrid float_t omega_fit(const std::size_t& s1, const std::size_t& s2, const float_t& T_star, const float_t& lambda_D, auto& omega_params) const
        {
            const float_t omega_nn = float_t(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(T_star*T_star))*std::log(omega_params(0, map_indices(s1, s2))*T_star*(float_t(1)-omega_params(1, map_indices(s1, s2))*std::exp(-1*omega_params(2, map_indices(s1, s2))*T_star))+float_t(1));
            return omega_nn;
        }

        // omega_11 collision integral (for neutrals)
        _sp_hybrid float_t omega_11(const std::size_t& s1, const std::size_t& s2, const float_t& T) const
        { return omega_fit(s1, s2, T, omega11_params); }

        // omega_11 collision integral (using coulomb shielding fit for ion-ion and ion-electron)
        _sp_hybrid float_t omega_11(const std::size_t& s1, const std::size_t& s2, const float_t& T_star, const float_t& lambda_D) const
        { return omega_fit(s1, s2, T_star, lambda_D, omega11_params); }
        
        // omega_22 collision integral (for neutrals)
        _sp_hybrid float_t omega_22(const std::size_t& s1, const std::size_t& s2, const float_t& T) const
        { return omega_fit(s1, s2, T, omega22_params); }

        // omega_22 collision integral (using coulomb shielding fit for ion-ion and ion-electron)
        _sp_hybrid float_t omega_22(const std::size_t& s1, const std::size_t& s2, const float_t& T_star, const float_t& lambda_D) const
        { return omega_fit(s1, s2, T_star, lambda_D, omega22_params); }

    };

    // NEED TO ADD THERMAL CONDUCTIVITY
    // ALSO CAN CLEAN UP BY BRANCHING OFF FOR DIFF T and FIT FUNCTION
    // ALSO FIX Ys (spade::ctrs::array<rtype, qface.nspecies()> Ys = fluid_state::get_Ys(qface);)
    // Gupta Viscous Model (Multicomponent Gas) <JRB | Implemented: 4-14-24 | Validated: TODO>
    template <typename dtype, typename gas_t> struct gupta_visc_t
    : public visc_law_interface_t<gupta_visc_t<dtype, gas_t>, omni::info_list_t<omni::info::value>>
    {
        typedef dtype float_t;

        using base_t = visc_law_interface_t<gupta_visc_t<dtype, gas_t>, omni::info_list_t<omni::info::value>>;
        using base_t::get_visc;
        using base_t::get_beta;
        using base_t::get_diffuse;
        using base_t::info_type;
        using fit_t = gupta_collision_fit_t<dtype, gas_t>;

        _sp_hybrid constexpr static std::size_t num_species(){return gas_t::nspecies();}

        /*
        // this function maps down a triangular matrix into a linear array
        _sp_hybrid constexpr static std::size_t map_indices(const std::size_t s1, const std::size_t s2)
        {
            return std::size_t(ns*std::min(s1, s2) + std::max(s2, s1) - (std::min(s1, s2)*(std::min(s1, s2)+1))/2);
        }*/

        // the fit parameter arrays, with n(n+1)/2 elements, corresponding to each unique combination of species
        //spade::ctrs::array<dtype, std::size_t((ns*(ns+1))/2)> A0, A1, A2, A3, B0, B1, B2, B3;

        const gas_t gas;
        const fit_t fit_fns;

        // constructor for gupta transport model
        gupta_visc_t(const gas_t& gas_in, const fit_t& fit_fns_in) : gas{gas_in}, fit_fns{fit_fns_in} { }
        
        // compute the viscosity
        _sp_hybrid float_t get_visc(const fluid_state::prim_chem_t<float_t, num_species()>& q) const
        {

            //std::cout << "TEST T: " << q.T() << "\n";

            // the output variable
            float_t mu_out = 0;

            // get the molar concentrations and density
            //const spade::ctrs::array<float_t, ns> Xr   = fluid_state::get_Xr(q, gas);
            const spade::ctrs::array<float_t, num_species()> Ys   = fluid_state::get_Ys(q);
            const spade::ctrs::array<float_t, num_species()> rhos = fluid_state::get_rhos(q, gas);

            // get the electron number density in (#/cm^3)
            float_t n_e = 0;
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
                const float_t mass_s = gas.mw_s[s1]/consts::Na_kmol;         

                // the denominator term
                float_t denominator = 0;
                
                // compute the collision terms at different controlling temperatures
                for (int s2 = 0; s2 < num_species(); s2++)
                {
                    // variables used to calulate transport properties
                    float_t T, T_star, lambda_D, omega_22, delta_22;

                    // if an electron is involved, use vib/elec temperature
                    if (gas.mw_s[s1] < 1 || gas.mw_s[s2] < 1)
                    {
                        T = q.Tv();
                        T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));
                        lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                    }
                    else
                    {
                        T = q.T();
                        T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));
                        lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                    }

                    // if interaction is neutral, use the neutral fit functions
                    if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                    {
                        omega_22 = fit_fns.omega_22(s1, s2, T);
                    }
                    // otherwise use the coulomb-shielding fit functions
                    else
                    {
                        omega_22 = fit_fns.omega_22(s1, s2, T_star, lambda_D);
                    }
                    delta_22 = float_t(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*T*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega_22;
                    //std::cout << "TEST omega22: " << omega_22 << "\n";
                    //std::cout << "TEST delta22: " << delta_22 << "\n";
                    /*

                    float_t delta2 = 0;
                    float_t omega2;

                    // an electron
                    if (gas.mw_s[s1] < 1)
                    {
                        if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                        {
                            // neutral-electron (vibrational/electronic controlling temperature)
                            omega2 = B3[map_indices(s1, s2)]*std::pow(q.Tv(), B0[map_indices(s1, s2)]*std::log(q.Tv())*std::log(q.Tv()) + B1[map_indices(s1, s2)]*std::log(q.Tv()) + B2[map_indices(s1, s2)]);
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
                            omega2 = B3[map_indices(s1, s2)]*std::pow(q.T(), B0[map_indices(s1, s2)]*std::log(q.T())*std::log(q.T()) + B1[map_indices(s1, s2)]*std::log(q.T()) + B2[map_indices(s1, s2)]);
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.T()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                        else if (gas.mw_s[s2] < 1)
                        {
                            // ion-electron (vibrational/electronic controlling temperature)
                            const dtype lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                            const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));

                            const dtype c2 = B0[map_indices(s1, s2)];
                            const dtype C2 = B1[map_indices(s1, s2)];
                            const dtype D2 = B2[map_indices(s1, s2)];

                            omega2 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.Tv()*q.Tv()))*std::log(D2*T_star*(dtype(1)-C2*std::exp(-1*c2*T_star))+dtype(1));
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.Tv()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                        else
                        {
                            // ion-ion (translational controlling temperature)
                            const dtype lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                            const dtype T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));

                            const dtype c2 = B0[map_indices(s1, s2)];
                            const dtype C2 = B1[map_indices(s1, s2)];
                            const dtype D2 = B2[map_indices(s1, s2)];
                            
                            omega2 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.T()*q.T()))*std::log(D2*T_star*(dtype(1)-C2*std::exp(-1*c2*T_star))+dtype(1));
                            delta2 = dtype(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*q.T()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega2;
                        }
                    } */

                    // add the result to the denominator term
                    denominator += Ys[s2]*gas.mw_si[s2]*delta_22;
                    //std::cout << "TEST: " << denominator << "\n";
                }

                // now add the result to the viscosity
                mu_out += mass_s*Ys[s1]*gas.mw_si[s1] / denominator;
            }

            return mu_out;
        }
        
        // compute the second viscosity
        _sp_hybrid dtype get_beta(const fluid_state::prim_chem_t<dtype, num_species()>& q) const
        {
            return -0.66666666667*this->get_visc(q);
        }
        
        // compute the species' diffusion coefficients
        _sp_hybrid spade::ctrs::array<dtype, num_species()> get_diffuse(const fluid_state::prim_chem_t<dtype, num_species()>& q) const
        {
            spade::ctrs::array<dtype, num_species()> diffuse_out;

            // get the molar concentrations and density
            //const spade::ctrs::array<dtype, ns> Xr = fluid_state::get_Xr(q, gas);
            const spade::ctrs::array<dtype, num_species()> Ys = fluid_state::get_Ys(q);
            const spade::ctrs::array<dtype, num_species()> rhos = fluid_state::get_rhos(q, gas);

            // compute the sum of molar concentrations
            dtype Ys_sum = 0;
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                Ys_sum +=  Ys[s1]*gas.mw_si[s1];
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

                        // variables used to calulate transport properties
                        float_t T, T_star, lambda_D, omega_11, delta_11, binary_diff;

                        // if an electron is involved, use vib/elec temperature
                        if (gas.mw_s[s1] < 1 || gas.mw_s[s2] < 1)
                        {
                            T = q.Tv();
                            T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));
                            lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                        }
                        else
                        {
                            T = q.T();
                            T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));
                            lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                        }

                        // if interaction is neutral, use the neutral fit functions
                        if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                        {
                            omega_11 = fit_fns.omega_11(s1, s2, T);
                        }
                        // otherwise use the coulomb-shielding fit functions
                        else
                        {
                            omega_11 = fit_fns.omega_11(s1, s2, T_star, lambda_D);
                        }
                        delta_11 = float_t(8)/3*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2] / (consts::pi*consts::Rgas_uni*T*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega_11;
                        binary_diff = consts::kSI*T/(q.p()*delta_11);

                        /*
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

                                const dtype c1 = A0[map_indices(s1, s2)];
                                const dtype C1 = A1[map_indices(s1, s2)];
                                const dtype D1 = A2[map_indices(s1, s2)];
                            
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

                                const dtype c1 = A0[map_indices(s1, s2)];
                                const dtype C1 = A1[map_indices(s1, s2)];
                                const dtype D1 = A2[map_indices(s1, s2)];
                            
                                omega1 = dtype(5)*std::pow(10, 15)*consts::pi*(lambda_D*lambda_D/(q.T()*q.T()))*std::log(D1*T_star*(dtype(1)-C1*std::exp(-1*c1*T_star))+dtype(1));
                            }
                            delta1 = dtype(8)/3*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/ (consts::pi*consts::Rgas_uni*q.T()*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10, -20)*omega1;
                            binary_diff = consts::kSI*q.T()/(q.p()*delta1);
                        } */

                        // add the term to the denominator
                        denominator += Ys[s2]*gas.mw_si[s2]/binary_diff;
                    }
                }

                // now compute and store the species diffusion coefficient
                diffuse_out[s1] = Ys_sum*Ys_sum*gas.mw_s[s1]*(1-Ys[s1])/denominator;
            }

            return diffuse_out;
        }

        // compute the translational & rotational thermal conductivity
        _sp_hybrid float_t get_kappa_tr(const fluid_state::prim_chem_t<float_t, num_species()>& q) const
        {
            // the translational and rotational conductivities
            float_t kappa_t, kappa_r = 0;

            // get the molar concentrations and density
            //const spade::ctrs::array<float_t, ns> Xr   = fluid_state::get_Xr(q, gas);
            const spade::ctrs::array<float_t, num_species()> Ys   = fluid_state::get_Ys(q);
            const spade::ctrs::array<float_t, num_species()> rhos = fluid_state::get_rhos(q, gas);

            // get the electron number density in (#/cm^3)
            float_t n_e = 0;
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
                // s1 cannot be an electron
                if (gas.mw_s[s1] >= 1)
                {
                    // species 1 molecular mass in kg
                    const float_t mass_s = gas.mw_s[s1]/consts::Na_kmol;         

                    // the denominator terms
                    float_t denominator_trans = 0;
                    float_t denominator_rot   = 0;
                    
                    // compute the collision terms at different controlling temperatures
                    for (int s2 = 0; s2 < num_species(); s2++)
                    {
                        // species 2 molecular mass in kg
                        const float_t mass_r = gas.mw_s[s2]/consts::Na_kmol; 
                        
                        // variables used to calulate transport properties
                        float_t T, T_star, lambda_D, omega_11, omega_22, delta_11, delta_22;

                        // if an electron is involved, use vib/elec temperature
                        if (gas.mw_s[s1] < 1 || gas.mw_s[s2] < 1)
                        {
                            T = q.Tv();
                            T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));
                            lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                        }
                        else
                        {
                            T = q.T();
                            T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));
                            lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                        }

                        // if interaction is neutral, use the neutral fit functions
                        if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                        {
                            omega_11 = fit_fns.omega_11(s1, s2, T);
                            omega_22 = fit_fns.omega_22(s1, s2, T);
                        }
                        // otherwise use the coulomb-shielding fit functions
                        else
                        {
                            omega_11 = fit_fns.omega_11(s1, s2, T_star, lambda_D);
                            omega_22 = fit_fns.omega_22(s1, s2, T_star, lambda_D);
                        }

                        delta_11 = float_t(8)/3*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*T*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega_11;
                        delta_22 = float_t(16)/5*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*T*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega_22;

                        // increment the denominator term for rotational thermal conductivity
                        denominator_rot += Ys[s2]*gas.mw_si[s2]*delta_11;

                        // increment the denominator term for translational thermal conductivity
                        if (gas.mw_s[s1] < 1 || gas.mw_s[s2] < 1)
                        {
                            denominator_trans += float_t(3.54)*Ys[s2]*gas.mw_si[s2]*delta_22;
                        }
                        else
                        {
                            const float_t a_sr = float_t(1)+(float_t(1)-mass_s/mass_r)*(float_t(0.45)-2.54*mass_s/mass_r)/((float_t(1)+mass_s/mass_r)*(float_t(1)+mass_s/mass_r));
                            denominator_trans += a_sr*Ys[s2]*gas.mw_si[s2]*delta_22;
                        }
                    }

                    // now add the result to the conductivities
                    kappa_t += float_t(15)/4*consts::kSI*Ys[s1]*gas.mw_si[s1] / denominator_trans;
                    kappa_r += consts::kSI*Ys[s1]*gas.get_cvr(s1) / (denominator_rot*consts::Na*float_t(1000));
                }
            }
            return kappa_t + kappa_r;
        }

        // compute the vibrational & electronic thermal conductivity
        _sp_hybrid float_t get_kappa_ve(const fluid_state::prim_chem_t<float_t, num_species()>& q) const
        {
            // the translational and rotational conductivities
            float_t kappa_ve = 0;

            // get the molar concentrations and density
            //const spade::ctrs::array<float_t, ns> Xr   = fluid_state::get_Xr(q, gas);
            const spade::ctrs::array<float_t, num_species()> Ys   = fluid_state::get_Ys(q);
            const spade::ctrs::array<float_t, num_species()> rhos = fluid_state::get_rhos(q, gas);

            // get the electron number density in (#/cm^3)
            float_t n_e = 0;
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                if (gas.mw_s[s1] < 1)
                {
                    n_e = std::pow(10, 6)*consts::Na_kmol*rhos[s1]/gas.mw_s[s1];
                    break;
                }
            }

            // now we want to compute the species' heat capacity at constant pressure
            spade::ctrs::array<float_t, num_species()> cp_ve;
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                cp_ve[s1] = gas.mw_s[s1] * (gas.get_cvv(s1, q.Tv()) + gas.get_cve(s1, q.Tv())) / (float_t(1e3) * consts::Na);
            }

            // now we can compute the vibrational/electronic thermal conductivity
            for (int s1 = 0; s1 < num_species(); s1++)
            {
                float_t denominator = 0;
                for (int s2 = 0; s2 < num_species(); s2++)
                {
                    // variables used to calulate transport properties
                    float_t T, T_star, lambda_D, omega_11, omega_22, delta_11, delta_22;

                    // if an electron is involved, use vib/elec temperature
                    if (gas.mw_s[s1] < 1 || gas.mw_s[s2] < 1)
                    {
                        T = q.Tv();
                        T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.Tv()));
                        lambda_D = sqrt(consts::kCGS*q.Tv()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                    }
                    else
                    {
                        T = q.T();
                        T_star = lambda_D/(consts::eCGS*consts::eCGS/(consts::kCGS*q.T()));
                        lambda_D = sqrt(consts::kCGS*q.T()/(4*consts::pi*n_e*consts::eCGS*consts::eCGS));
                    }

                    // if interaction is neutral, use the neutral fit functions
                    if (gas.charge_s[s1] == 0 || gas.charge_s[s2] == 0)
                    {
                        omega_11 = fit_fns.omega_11(s1, s2, T);
                    }
                    // otherwise use the coulomb-shielding fit functions
                    else
                    {
                        omega_11 = fit_fns.omega_11(s1, s2, T_star, lambda_D);
                    }

                    delta_11 = float_t(8)/3*sqrt(2*gas.mw_s[s1]*gas.mw_s[s2]/(consts::pi*consts::Rgas_uni*T*(gas.mw_s[s1]+gas.mw_s[s2])))*std::pow(10,-20)*omega_11;
                    denominator += Ys[s2]*gas.mw_si[s2]*delta_11;
                }

                kappa_ve += Ys[s1]*gas.mw_si[s1]*cp_ve[s1] / denominator;
            }

            return kappa_ve;
        }

    };

    // Initialization Function for the Collision Integral Fit Parameters <JRB | Implemented: 4-14-24 | Validated: TODO>
    /*
	template<typename dtype, const std::size_t ns, fluid_state::is_multicomponent_gas_type gas_t>
	static void import_gupta_collision_integral_data(const std::string& fname, const std::vector<std::string>& speciesNames, gas_t& gas, gupta_visc_t<dtype, ns, gas_t>& gupta_visc)
	{
		// Open the input file
        std::ifstream infile;
		std::string full_fname = std::getenv("SPADE") + std::string("/src/navier-stokes/speciesCollision/") + fname;
        infile.open(full_fname);
		
		if (infile)
		{
			// variables to temporarily store data
			std::string species1, species2;
            dtype A0, A1, A2, A3, B0, B1, B2, B3;
            spade::ctrs::array<std::size_t, ns*ns> unique_indices, unique_indices_check;

			// sweep entire file
			while (true)
			{
                // collision file is formatted: collision species, omega11 params (A), omega22 params (B)
				infile >> species1 >> species2 >> A0 >> A1 >> A2 >> A3 >> B0 >> B1 >> B2 >> B3;
				
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
                                gupta_visc.B0[gupta_visc.map_indices(s1, s2)] = B0;
                                gupta_visc.B1[gupta_visc.map_indices(s1, s2)] = B1;
                                gupta_visc.B2[gupta_visc.map_indices(s1, s2)] = B2;
                                gupta_visc.B3[gupta_visc.map_indices(s1, s2)] = B3;

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

            // the coulomb shieliding collision fit parameters
            constexpr static dtype c1_att = 0.0313, C1_att = -0.476, D1_att = 0.784;
            constexpr static dtype c1_rep = 0.0106, C1_rep = 0.138,  D1_rep = 0.765;
            constexpr static dtype c2_att = 0.0377, C2_att = -0.146, D2_att = 1.262;
            constexpr static dtype c2_rep = 0.0274, C2_rep = 0.157,  D2_rep = 1.235;

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
                    else
                    {
                        // if the collision is not neutral, use the coulomb sheilding fit parameters
                        if (spade::utils::sign(gas.charge_s[s1]) != spade::utils::sign(gas.charge_s[s2])) 
                        { 
                            A0 = c1_att; A1 = C1_att; A2 = D1_att; A3 = dtype(1);
                            B0 = c2_att; B1 = C2_att; B2 = D2_att; B3 = dtype(1);
                        }
                        else 
                        { 
                            A0 = c1_rep; A1 = C1_rep; A2 = D1_rep; A3 = dtype(1);
                            B0 = c2_rep; B1 = C2_rep; B2 = D2_rep; B3 = dtype(1);
                        }
                        gupta_visc.A0[gupta_visc.map_indices(s1, s2)] = A0;
                        gupta_visc.A1[gupta_visc.map_indices(s1, s2)] = A1;
                        gupta_visc.A2[gupta_visc.map_indices(s1, s2)] = A2;
                        gupta_visc.A3[gupta_visc.map_indices(s1, s2)] = A3;
                        gupta_visc.B0[gupta_visc.map_indices(s1, s2)] = B0;
                        gupta_visc.B1[gupta_visc.map_indices(s1, s2)] = B1;
                        gupta_visc.B2[gupta_visc.map_indices(s1, s2)] = B2;
                        gupta_visc.B3[gupta_visc.map_indices(s1, s2)] = B3;
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

    */ // TO REMOVE
}