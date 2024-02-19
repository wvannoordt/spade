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
        t.get_diffuse();  // cp * mu/Pr
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
        : visc{visc_in}, beta{-2.0*visc_in/3.0}, prandtl_inv{dtype(1.0)/prandtl_in}
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

    template <typename laminar_t, typename turb_t>
    struct sgs_visc_t
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
                lam_all.beta  - 0.66666666667*mu_t,
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
}