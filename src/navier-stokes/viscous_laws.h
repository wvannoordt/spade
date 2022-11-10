#pragma once

#include "navier-stokes/fluid_state.h"
#include <concepts>


namespace spade::viscous_laws
{
    template <class T> concept state_independent_viscosity = std::floating_point<typename T::value_type> && requires(T t)
    {
        typename T::value_type;
        t.get_visc();
        t.get_beta();
        t.get_conductivity();
    };
    
    template <class T> concept state_dependent_viscosity = std::floating_point<typename T::value_type>
    && requires(T t, const fluid_state::prim_t<typename T::value_type>& q)
    {
        typename T::value_type;
        t.get_visc(q);
        t.get_beta(q);
        t.get_conductivity(q);
    };
    
    template <class T> concept viscous_law = state_dependent_viscosity<T> || state_independent_viscosity<T>;
    
    template <typename dtype> struct constant_viscosity_t
    {
        typedef dtype value_type;
        constant_viscosity_t(const dtype& visc_in, const dtype& prandtl_in) { this->visc = visc_in; this->beta = -0.66666666667*this->visc; prandtl = prandtl_in;}
        constant_viscosity_t(void) { this->visc = dtype(); this->beta = -0.66666666667*this->visc; }
        
        template <fluid_state::is_state_type state_t> dtype get_visc(const state_t& q) const
        {
            return visc;
        }
        
        template <fluid_state::is_state_type state_t> dtype get_beta(const state_t& q) const
        {
            return beta;
        }
        
        template <fluid_state::is_state_type state_t> dtype get_conductivity(const state_t& q) const
        {
            return visc/prandtl;
        }
        
        dtype get_visc(void) const
        {
            return visc;
        }
        
        dtype get_beta(void) const
        {
            return beta;
        }
        
        dtype get_conductivity(void) const
        {
            return visc/prandtl;
        }
        
        dtype visc;
        dtype prandtl;
        dtype beta;
    };
    
    template <typename dtype> struct power_law_t
    {
            typedef dtype value_type;
            
            dtype mu_ref;
            dtype T_ref;
            dtype power;
            dtype prandtl;
            
            power_law_t(const dtype& mu_ref_in, const dtype& T_ref_in, const dtype& power_in, const dtype& prandtl_in)
            {
                mu_ref  = mu_ref_in;
                T_ref   = T_ref_in;
                power   = power_in;
                prandtl = prandtl_in;
            }
            
            template <fluid_state::is_state_type state_t> dtype get_visc(const state_t& q) const
            {
                return mu_ref*std::pow(q.T()/T_ref, power);
            }
            
            template <fluid_state::is_state_type state_t> dtype get_beta(const state_t& q) const
            {
                return -0.66666666667*this->get_visc(q);
            }
            
            template <fluid_state::is_state_type state_t> dtype get_conductivity(const state_t& q) const
            {
                return this->get_visc(q)/prandtl;
            }
    };
}