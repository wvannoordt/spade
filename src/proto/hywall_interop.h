#pragma once

#include <concepts>
#include <cstring>
#include <iostream>
#include <vector>

#include "HyWall.h"
#include "PTL.h"

#include "navier-stokes/fluid_state.h"

#include "proto/hywall_read_inputs.h"

namespace spade::proto
{
    template <typename array_t, typename gas_t, typename real_t = double>
    requires (fluid_state::is_state_type<typename array_t::alias_type> && fluid_state::state_convertible<typename array_t::alias_type, fluid_state::prim_t<real_t>, gas_t>)
    struct hywall_binding_t
    {
        std::size_t num_points;
        const gas_t* gas;
        std::vector<real_t> in_prims;
        
        std::vector<real_t> in_dist;
        std::vector<real_t> in_x;
        std::vector<real_t> in_y;
        std::vector<real_t> in_z;
        std::vector<real_t> in_rho;
        std::vector<real_t> in_mu;
        std::vector<real_t> in_momrhs;
        std::vector<real_t> in_dpdx;
        
        std::vector<real_t> aux_strainrate;
        std::vector<real_t> aux_sensorpremult;
        
        std::vector<real_t> out_vort;
        std::vector<real_t> out_tau;
        std::vector<real_t> out_qw;
        std::vector<real_t> out_fail;
        
        using state_t = typename array_t::alias_type;
        
        hywall_binding_t(const array_t& prim_in, const gas_t& gas_in)
        {
            gas = &gas_in;
        }
        
        void read(PTL::PropertySection& input_section)
        {
            hywall_read_inputs(input_section, HyWall::settings);
        }
        
        void init(const array_t& q)
        {
            //note: storage in in_prims is p p p p p p p p p p p p u u u u u u u u u u u u u v v v v v v v...
            std::size_t num_wm_points = 100;
            this->set_size(num_wm_points);
            HyWall::SetDomainSize(num_wm_points);
            HyWall::DefineVariables();
            HyWall::PassFlowfieldVariables(&in_prims[0], num_wm_points);
            HyWall::PassVariable("in:distance",    &in_dist  [0]);
            HyWall::PassVariable("in:x",           &in_x     [0]);
            HyWall::PassVariable("in:y",           &in_y     [0]);
            HyWall::PassVariable("in:z",           &in_z     [0]);
            HyWall::PassVariable("in:rho",         &in_rho   [0]);
            HyWall::PassVariable("in:mu_lam",      &in_mu    [0]);
            HyWall::PassVariable("in:momRHS",      &in_momrhs[0]);
            HyWall::PassVariable("in:dpdx",        &in_dpdx  [0]);
        	if (HyWall::settings.enableTransitionSensor)
        	{
        		HyWall::PassVariable("aux:strain_rate",    &aux_strainrate   [0]);
        		HyWall::PassVariable("aux:sensor_preMult", &aux_sensorpremult[0]);
        	}
        	
        	HyWall::PassVariable("out:vorticity",    &out_vort[0]);
        	HyWall::PassVariable("out:tau",          &out_tau [0]);
        	HyWall::PassVariable("out:heatflux",     &out_qw  [0]);
        	HyWall::PassVariable("out:failurelevel", &out_fail[0]);
            HyWall::Allocate();
        }
        
        void set_size(const std::size_t size_in)
        {
            num_points = size_in;
            in_prims.resize(num_points*6);
            in_dist.resize(num_points);
            in_x.resize(num_points);
            in_y.resize(num_points);
            in_z.resize(num_points);
            in_rho.resize(num_points);
            in_mu.resize(num_points);
            in_momrhs.resize(num_points);
            in_dpdx.resize(num_points);
            aux_strainrate.resize(num_points);
            aux_sensorpremult.resize(num_points);
            out_vort.resize(num_points);
            out_tau.resize(num_points);
            out_qw.resize(num_points);
            out_fail.resize(num_points);
        }
        
        std::size_t size() const {return num_points;}
        
        void set_dt(const real_t& dt)
        {
            HyWall::SetTimeStep(dt);
        }
        
        ~hywall_binding_t()
        {
            HyWall::Finalize();
        }
    };
}