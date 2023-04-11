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

    template <typename derived_t, omni::info::is_omni_info... infos_t>
    struct visc_law_interface_t
    {
        derived_t&       self()       {return *static_cast<      derived_t*>(this);}
        const derived_t& self() const {return *static_cast<const derived_t*>(this);}

        using info_type = omni::info_list_t<infos_t...>;

        auto get_visc        (const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_visc        (args...);}, input);
        }

        auto get_beta        (const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_beta        (args...);}, input);
        }

        auto get_diffuse(const auto& input) const
        {
            return omni::invoke_call(info_type(), [&](const auto&... args){return this->self().get_diffuse(args...);}, input);
        }
    };
    
    template <typename dtype> struct constant_viscosity_t
    : public visc_law_interface_t<constant_viscosity_t<dtype>>
    {
        using base_t = visc_law_interface_t<constant_viscosity_t<dtype>>;
        using base_t::get_visc;
        using base_t::get_beta;
        using base_t::get_diffuse;
        using base_t::info_type;

        typedef dtype value_type;
        constant_viscosity_t(const dtype& visc_in, const dtype& prandtl_in)
        : visc{visc_in}, beta{-2.0*visc_in/3.0}, prandtl{prandtl_in}
        { }
        
        dtype get_visc() const
        {
            return visc;
        }
        
        dtype get_beta() const
        {
            return beta;
        }
        
        dtype get_diffuse() const
        {
            return visc/prandtl;
        }
        
        dtype visc;
        dtype prandtl;
        dtype beta;
    };
    
    template <typename dtype> struct power_law_t
    : public visc_law_interface_t<power_law_t<dtype>, omni::info::value>
    {
            typedef dtype value_type;
            dtype mu_ref;
            dtype T_ref;
            dtype power;
            dtype prandtl;

            using base_t = visc_law_interface_t<power_law_t<dtype>, omni::info::value>;
            using base_t::get_visc;
            using base_t::get_beta;
            using base_t::get_diffuse;
            using base_t::info_type;
            
            power_law_t(const dtype& mu_ref_in, const dtype& T_ref_in, const dtype& power_in, const dtype& prandtl_in)
            : mu_ref{mu_ref_in}, T_ref{T_ref_in}, power{power_in}, prandtl{prandtl_in}
            { }
            
            template <fluid_state::is_state_type state_t> dtype get_visc(const state_t& q) const
            {
                return mu_ref*std::pow(q.T()/T_ref, power);
            }
            
            template <fluid_state::is_state_type state_t> dtype get_beta(const state_t& q) const
            {
                return -0.66666666667*this->get_visc(q);
            }
            
            template <fluid_state::is_state_type state_t> dtype get_diffuse(const state_t& q) const
            {
                return this->get_visc(q)/prandtl;
            }
    };
    
    template <typename state_t, typename visc_func_t, typename beta_func_t, typename cond_func_t> struct udf_t
    : public visc_law_interface_t<udf_t<state_t, visc_func_t, beta_func_t, cond_func_t>, omni::info::value>
    {
            typedef typename state_t::value_type value_type;

            using base_t = visc_law_interface_t<udf_t<state_t, visc_func_t, beta_func_t, cond_func_t>, omni::info::value>;
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
            
            value_type get_visc(const state_t& q) const
            {
                return vfunc(q);
            }
            
            value_type get_beta(const state_t& q) const
            {
                return bfunc(q);
            }
            
            value_type get_diffuse(const state_t& q) const
            {
                return cfunc(q);
            }
    };
}