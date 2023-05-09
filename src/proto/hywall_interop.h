#pragma once

#include <cstring>
#include <iostream>
#include <vector>

#include "HyWall.h"

#include "core/config.h"
#include "navier-stokes/fluid_state.h"
#include "core/sdf_binding.h"

#include "proto/hywall_read_inputs.h"

namespace spade::proto
{
    template <typename array_t, typename rhs_t, typename gas_t, typename real_t = double>
    requires (fluid_state::is_state_type<typename array_t::alias_type> && fluid_state::state_convertible<typename array_t::alias_type, fluid_state::prim_t<real_t>, gas_t>)
    struct hywall_binding_t
    {
        using array_state_t = typename array_t::alias_type;
        using coord_type    = typename array_t::grid_type::coord_type;
        std::vector<grid::face_idx_t> wm_faces;
        std::vector<grid::cell_idx_t> wm_samples;
        std::vector<ctrs::v3d<coord_type>> wm_nvec_comps;
        std::vector<ctrs::v3d<coord_type>> wm_tvec;
        std::size_t num_points;
        int sampling_distance;
        int nx, nz;
        real_t dx;
        const gas_t* gas;
        const array_t::grid_type::group_type& group;
        
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
        std::vector<real_t> out_tau2;
        std::vector<real_t> out_qw2;
        std::vector<real_t> out_fail;
        
        using state_t = typename array_t::alias_type;
        
	  hywall_binding_t(const array_t& prim_in, const rhs_t& rhs, const gas_t& gas_in, const int& sampling_distance_in)
        : group{prim_in.get_grid().group()}, sampling_distance{sampling_distance_in}
        {
            gas = &gas_in;
            HyWall::Initialize(prim_in.get_grid().group().get_channel(), 0);
            dx = prim_in.get_grid().get_dx(1);
            nx = prim_in.get_grid().get_num_cells(0);
            nz = prim_in.get_grid().get_num_cells(2);
        }
        
        void read(auto& input_section)
        {
            hywall_read_inputs(input_section, HyWall::settings);
        }
        
        std::size_t enumerate_points(const array_t& q, const bound_box_t<bool, array_t::grid_type::dim()>& walls)
        {
            const std::size_t dim = array_t::grid_type::dim();
            const auto& grid = q.get_grid();
            const std::size_t nlb_loc = grid.get_num_local_blocks();
            for (auto lb_loc: range(0, nlb_loc))
            {
                std::size_t lb_glob = grid.get_partition().get_global_block(lb_loc);
                const auto& idomain = grid.is_domain_boundary(lb_glob);
                for (auto bndy: range(0,2)*range(0,dim))
                {
                    int pm  = bndy[0];
                    int xyz = bndy[1];
                    ctrs::v3d<coord_type> nvec_comp = 0;
                    nvec_comp[xyz] = 1-2*pm;
                    if (idomain(xyz, pm) && walls(xyz, pm))
                    {
                        bound_box_t<int, 3> cell_box;
                        for (auto d: range(0, 3))
                        {
                            cell_box.min(d) = 0;
                            cell_box.max(d) = grid.get_num_cells(d);
                        }
                        if (pm == 0) cell_box.max(xyz) = 1;
                        if (pm == 1) cell_box.min(xyz) = grid.get_num_cells(xyz)-1;
                        auto rg0 = range(cell_box.min(0), cell_box.max(0));
                        auto rg1 = range(cell_box.min(1), cell_box.max(1));
                        auto rg2 = range(cell_box.min(2), cell_box.max(2));
                        for (auto i: rg0*rg1*rg2)
                        {
                            const int sampl_dist = sampling_distance;
                            grid::cell_idx_t ijk(i[0], i[1], i[2], lb_loc);
                            grid::face_idx_t ijkf = grid::cell_to_face(ijk, xyz, pm);
                            wm_faces.push_back(ijkf);
                            ijk[xyz] -= (2*pm - 1)*sampl_dist;
                            wm_samples.push_back(ijk);
                            wm_nvec_comps.push_back(nvec_comp);
                        }
                    }
                }
            }
            this->set_size(wm_faces.size());
            for (auto n: range(0, wm_faces.size()))
            {
                auto ijk   = wm_samples[n];
                auto ijkF  = wm_faces[n];
                
                auto xc = grid.get_coords(ijk);
                auto xf = grid.get_coords(ijkF);
                auto dx = xc;
                dx -= xf;
                auto dist = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
                in_dist[n]     = dist;
                in_x[n]        = xf[0];
                in_y[n]        = xf[1];
                in_z[n]        = xf[2];
            }
            out_tau2 = out_tau;
            out_qw2  = out_qw;
            return wm_faces.size();
        }
        
        void init(const array_t& q, const bound_box_t<bool, array_t::grid_type::dim()>& walls)
        {
            //note: storage in in_prims is p p p p p p p p p p p p u u u u u u u u u u u u u v v v v v v v...
            std::size_t num_wm_points = this->enumerate_points(q, walls);
            if (this->size() == 0) return;
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
        
        void set_size(const std::size_t& size_in)
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
            wm_tvec.resize(num_points);
        }

        real_t avg(const auto& vec) const
        {
            std::size_t num = group.sum(vec.size());
            real_t sum = 0.0;
            for (const auto& p: vec) sum += p;
            sum = group.sum(sum);
            return sum/num;
        }

        real_t mean_tau() const
        {
            return this->avg(out_tau);
        }

        real_t mean_qw () const
        {
            return this->avg(out_qw);
        }
        
        std::size_t size() const {return num_points;}
        
        void set_dt(const real_t& dt)
        {
            HyWall::SetTimeStep(dt);
        }
        
        template <typename visc_t> void sample(const array_t& q, const visc_t& visc)
        {
            if (this->size() == 0) return;
            for (auto n: range(0, wm_faces.size()))
            {
                auto ijk = wm_samples[n];
                
                array_state_t state;
                for (auto j: range(0, state.size())) state[j] = q(j, ijk[0], ijk[1], ijk[2], ijk[3]);
                fluid_state::prim_t<real_t> prim;
                fluid_state::cons_t<real_t> cons;
                fluid_state::convert_state(state, prim, *gas);
                fluid_state::convert_state(state, cons, *gas);
                
                ctrs::v3d<real_t> tvec(prim.u(), 0.0, prim.w());
                tvec /= ctrs::array_norm(tvec);
                wm_tvec[n] = tvec;
                
                real_t p = prim.p();
                real_t T = prim.T();
                real_t u = std::sqrt(prim.u()*prim.u() + prim.w()*prim.w());
                real_t v = 0.0;
                real_t w = 0.0;
                real_t nu = 0.0;
                in_prims[0*wm_faces.size() + n]  = p;
                in_prims[1*wm_faces.size() + n]  = u;
                in_prims[2*wm_faces.size() + n]  = v;
                in_prims[3*wm_faces.size() + n]  = w;
                in_prims[4*wm_faces.size() + n]  = T;
                in_prims[5*wm_faces.size() + n]  = nu;
                
                in_mu [n]    = visc.get_visc(prim);
                in_momrhs[n] = 0.0;
                in_dpdx[n]   = 0.0;
                in_rho[n]    = cons.rho();
                
                aux_strainrate[n] = 1.0;
                aux_sensorpremult[n] = 1.0;
            }
        }
        
        void solve()
        {
            if (this->size() == 0) return;
            HyWall::Solve();
        }

        void smooth()
        {
            out_tau2 = out_tau;
            out_qw2  = out_qw;
            std::size_t idx = 0;
            for (auto& i_f: wm_faces)
            {
                if ((i_f.i() > 0) && (i_f.i() < nx-1) && (i_f.k() > 0) && (i_f.k() < nz - 1))
                {
                    auto s0 = idx;
                    auto s1 = idx + 1;
                    auto s2 = idx - 1;
                    auto s3 = idx + nx;
                    auto s4 = idx - nz;
                    auto tauval = 0.2*(out_tau2[s0] + out_tau2[s1] + out_tau2[s2] + out_tau2[s3] + out_tau2[s4]);
                    auto  qwval = 0.2*(out_qw2 [s0] + out_qw2 [s1] + out_qw2 [s2] + out_qw2 [s3] + out_qw2 [s4]);
                    out_tau[idx] = tauval;
                    out_qw [idx] = qwval;
                }
                idx++;
            }
        }
        
        void apply_flux(rhs_t& rhs) const
        {
            if (this->size() == 0) return;
            for (auto n: range(0, wm_faces.size()))
            {
                auto ijkF = wm_faces[n];
                auto ijk_l = grid::face_to_cell(ijkF, 0);
                auto ijk_r = grid::face_to_cell(ijkF, 1);

                fluid_state::flux_t<real_t> flux;
                flux.continuity() = 0.0;
                flux.energy()     = -out_qw[n] *wm_nvec_comps[n][1];
                flux.x_momentum() =  out_tau[n]*wm_tvec[n][0]*wm_nvec_comps[n][1];
                flux.y_momentum() =  0.0;
                flux.z_momentum() =  out_tau[n]*wm_tvec[n][2]*wm_nvec_comps[n][1];

                //This is correct 2023-03-30
                for (int iii = 0; iii < flux.size(); ++iii)
                {
                    rhs(iii, ijk_l) += flux[iii]/dx;
                    rhs(iii, ijk_r) -= flux[iii]/dx;
                }
            }
        }
        
        ~hywall_binding_t()
        {
            HyWall::Finalize();
        }
    };
}
