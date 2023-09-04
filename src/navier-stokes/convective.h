#pragma once

#include "core/config.h"
#include "core/ctrs.h"
#include "omni/omni.h"

#include "core/finite_diff.h"

#include "navier-stokes/fluid_state.h"

namespace spade::convective
{       
    template <typename gas_t> struct totani_lr
    {
        using float_t       = typename gas_t::value_type;
        using output_type   = fluid_state::flux_t<float_t>;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using sub_info_type = typename gas_t::info_type;
        using info_type     = omni::info_union<own_info_type, sub_info_type>;
        using omni_type     = omni::prefab::lr_t<info_type>;
        
        const gas_t& gas;
        
        totani_lr(const gas_t& gas_in) : gas{gas_in}{}
        
        
        output_type operator() (const auto& input_data) const
        {
            output_type output;
            const auto& ql       = omni::access<omni::info::value> (input_data.cell(0_c));
            const auto& qr       = omni::access<omni::info::value> (input_data.cell(1_c));
            const auto& normal_l = omni::access<omni::info::metric>(input_data.cell(0_c));
            const auto& normal_r = omni::access<omni::info::metric>(input_data.cell(1_c));
            
            float_t un_l = normal_l[0]*ql.u()+normal_l[1]*ql.v()+normal_l[2]*ql.w();
            float_t un_r = normal_r[0]*qr.u()+normal_r[1]*qr.v()+normal_r[2]*qr.w();
            
            float_t rho_l = ql.p()/(gas.get_R(input_data.cell(0_c))*ql.T());
            float_t rho_r = qr.p()/(gas.get_R(input_data.cell(1_c))*qr.T());
            
            float_t e_l = ql.p()/(rho_l*(gas.get_gamma(input_data.cell(0_c))-float_t(1.0)));
            float_t e_r = qr.p()/(rho_r*(gas.get_gamma(input_data.cell(1_c))-float_t(1.0)));
            
            float_t c = float_t(0.25)*(rho_l+rho_r)*(un_l+un_r);
            output.continuity() = c;
            output.energy()     = float_t(0.5)*c*(e_l + e_r + ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w()) + float_t(0.5)*(un_l*qr.p() + un_r*ql.p());
            output.x_momentum() = float_t(0.5)*c*(ql.u()+qr.u()) + float_t(0.5)*(normal_l[0]*ql.p()+normal_r[0]*qr.p());
            output.y_momentum() = float_t(0.5)*c*(ql.v()+qr.v()) + float_t(0.5)*(normal_l[1]*ql.p()+normal_r[1]*qr.p());
            output.z_momentum() = float_t(0.5)*c*(ql.w()+qr.w()) + float_t(0.5)*(normal_l[2]*ql.p()+normal_r[2]*qr.p());
            return output;
        }
    };
    
    template <typename gas_t, const int order>
    struct cent_keep_scheme_t
    {
        static_assert(order <= 8, "only accurate up to 8th order because of precision issues");
        constexpr static int half_wid = order/2;
        
        using float_t       = typename gas_t::value_type;
        using output_type   = fluid_state::flux_t<float_t>;
        using own_info_type = omni::info_list_t<omni::info::value, omni::info::metric>;
        using sub_info_type = typename gas_t::info_type;
        using info_type     = omni::info_union<own_info_type, sub_info_type>;
        using x0_t          = omni::offset_t< -(2*half_wid-1), 0, 0>;
        using dx_t          = omni::offset_t<               2, 0, 0>;
        using cells_t       = omni::pattern_t<grid::face_centered, x0_t, dx_t, order, info_type>;
        using faces_t       = omni::prefab::face_mono_t<omni::info_list_t<omni::info::metric>>;
        using omni_type     = omni::stencil_union<cells_t, faces_t>;
        
        const gas_t gas;
        
        cent_keep_scheme_t(const gas_t& gas_in) : gas{gas_in}{}
        
        _sp_hybrid output_type operator() (const auto& input_data) const
        {
            output_type output;
            const auto nface = omni::access<omni::info::metric>(input_data.face(0_c));
            ctrs::array<float_t, order> u_normal, rho, eng;
            algs::static_for<0, order>([&](const auto& ii)
            {
                udci::idx_const_t<ii.value> idx;
                const auto& c_cell = input_data.cell(idx);
                const auto& q = omni::access<omni::info::value >(c_cell);
                const auto& n = omni::access<omni::info::metric>(c_cell);
                u_normal[idx] = q.u()*n[0] + q.v()*n[1] + q.w()*n[2];
                rho[idx]      = q.p()/(gas.get_R(c_cell)*q.T());
                eng[idx]      = q.p()/(rho[idx]*(gas.get_gamma(c_cell)-float_t(1.0)));
            });
            float_t c   = float_t(0.0); // mass conservation
            float_t mx  = float_t(0.0); // x-momentum (advection)
            float_t my  = float_t(0.0); // y-momentum (advection)
            float_t mz  = float_t(0.0); // z-momentum (advection)
            float_t g   = float_t(0.0); // pressure gradient
            float_t k   = float_t(0.0); // kinetic energy flux
            float_t i   = float_t(0.0); // internal energy flux
            float_t p   = float_t(0.0); // pressure diffusion term
            
            // functors for stencil point calculations
            const auto ar_access = [&](const auto& arr){ return [&](const auto& iii) { return arr[iii.value]; }; };
            const auto st_access = [&](const auto& lam)
            {
                return [&](const auto& iii)
                {
                    const auto& st = omni::access<omni::info::value>(input_data.cell(udci::idx_const_t<iii.value>()));
                    return lam(st);
                };
            };
            algs::static_for<1, half_wid+1>([&](const auto& ii)
            {
                constexpr auto a_coeff = finite_diff::centered_finite_diff_t<ii.value, order, float_t>::value;
                algs::static_for<0, ii.value>([&](const auto& jj)
                {
                    udci::idx_const_t<half_wid - 1 - jj.value> i0;
                    udci::idx_const_t<i0.value + ii.value> i1;
                    const auto& q0 = omni::access<omni::info::value>(input_data.cell(udci::idx_const_t<i0.value>()));
                    const auto& q1 = omni::access<omni::info::value>(input_data.cell(udci::idx_const_t<i1.value>()));
                    //local averaging operators
                    const auto avg0 = [&](const auto& func)                     { return float_t(0.5)*(func(i0) + func(i1)); };
                    const auto avg1 = [&](const auto& func0, const auto& func1) { return float_t(0.25)*(func0(i0) + func0(i1))*(func1(i1) + func1(i0)); };
                    const auto avg2 = [&](const auto& func0, const auto& func1) { return float_t(0.5)*(func0(i0)*func1(i1) + func0(i1)*func1(i0)); };
                    
                    const auto c_loc = a_coeff*avg1(ar_access(rho), ar_access(u_normal));
                    c  +=   c_loc;
                    mx +=   c_loc*avg0(st_access([](const auto& q){return q.u();}));
                    my +=   c_loc*avg0(st_access([](const auto& q){return q.v();}));
                    mz +=   c_loc*avg0(st_access([](const auto& q){return q.w();}));
                    g  += a_coeff*avg0(st_access([](const auto& q){return q.p();}));
                    k  += float_t(0.5)*c_loc*(q0.u()*q1.u() + q0.v()*q1.v() + q0.w()*q1.w());
                    i  += c_loc*avg0(ar_access(eng));
                    p  += a_coeff*avg2(ar_access(u_normal), st_access([](const auto& q){return q.p();}));
                });
            });
            
            output.continuity() = float_t(2.0)*c;
            output.energy()     = float_t(2.0)*(k + i + p);
            output.x_momentum() = float_t(2.0)*(mx + g*nface[0]);
            output.y_momentum() = float_t(2.0)*(my + g*nface[1]);
            output.z_momentum() = float_t(2.0)*(mz + g*nface[2]);
            return output;
        }
    };
    
    template <const int order, typename gas_t>
    auto cent_keep(const gas_t& gas)
    {
        static_assert(order == 2*(order/2), "nominal order for centered scheme must be even");
        return cent_keep_scheme_t<gas_t, order>(gas);
    }
    
    template <typename gas_t, typename scale_t> struct pressure_diss_lr
    {
        typedef typename gas_t::value_type dtype;
        typedef fluid_state::flux_t<dtype> output_type;
        using sub_omni_type = omni::prefab::mono_t<grid::face_centered, typename scale_t::info_type>;
        using sup_omni_type = omni::prefab::lr_t<omni::info_list_t<omni::info::value, omni::info::metric>>;
        using omni_type     = omni::stencil_union<sub_omni_type, sup_omni_type>;

        const gas_t   gas;
        const scale_t scale;
        dtype eps_p, eps_T;
        
        pressure_diss_lr(const gas_t& gas_in, const scale_t& scale_in, const dtype& eps_p_in)
        : gas{gas_in}, scale{scale_in}, eps_p{eps_p_in}, eps_T{eps_p_in} {}

        pressure_diss_lr(const gas_t& gas_in, const scale_t& scale_in, const dtype& eps_p_in, const dtype& eps_T_in)
        : gas{gas_in}, scale{scale_in}, eps_p{eps_p_in}, eps_T{eps_T_in} {}
        
        _sp_hybrid output_type operator() (const auto& input) const
        {
            output_type output;
            const auto& ql       = omni::access<omni::info::value >(input.cell(0_c));
            const auto& qr       = omni::access<omni::info::value >(input.cell(1_c));
            const auto& normal_l = omni::access<omni::info::metric>(input.cell(0_c));
            const auto& normal_r = omni::access<omni::info::metric>(input.cell(1_c));
            
            const dtype delta_p = ql.p()-qr.p();
            const dtype delta_T = ql.T()-qr.T();
            const dtype pavg = 0.5*(ql.p() + qr.p());
            const dtype Tavg = 0.5*(ql.T() + qr.T());
            const dtype inv_RT = 1.0/(0.5*(gas.get_R(input.cell(0_c))+gas.get_R(input.cell(1_c)))*Tavg);
            const dtype p_inv_rt2 = pavg*inv_RT/Tavg;
            const dtype gamma =  0.5*(gas.get_gamma(input.cell(0_c))+gas.get_gamma(input.cell(1_c)));
            const dtype k = 0.5*(ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w());
            const dtype u = 0.5*(ql.u()+qr.u());
            const dtype v = 0.5*(ql.v()+qr.v());
            const dtype w = 0.5*(ql.w()+qr.w());
            const auto  scl      = scale.get_sensor(input.face(0_c));
            output.continuity()  = scl*eps_p*delta_p*inv_RT;
            output.energy()      = scl*eps_p*delta_p*(k*inv_RT + 1.0/gamma);
            output.x_momentum()  = scl*eps_p*delta_p*u*inv_RT;
            output.y_momentum()  = scl*eps_p*delta_p*v*inv_RT;
            output.z_momentum()  = scl*eps_p*delta_p*w*inv_RT;
            
            output.continuity() -= scl*eps_T*delta_T*p_inv_rt2;
            output.energy()     -= scl*eps_T*delta_T*p_inv_rt2*(0.5*k);
            output.x_momentum() -= scl*eps_T*delta_T*p_inv_rt2*u;
            output.y_momentum() -= scl*eps_T*delta_T*p_inv_rt2*v;
            output.z_momentum() -= scl*eps_T*delta_T*p_inv_rt2*w;
            
            return output;
        }
    };

    //todo: generalize to other weno schemes
    enum weno_smooth_indicator
    {
        enable_smooth,
        disable_smooth
    };

    template <typename flux_func_t, const weno_smooth_indicator use_smooth = enable_smooth>
    struct weno_t
    {
        using float_t       = typename flux_func_t::float_t;
        using output_type   = fluid_state::flux_t<float_t>;
        using info_type     = typename flux_func_t::info_type;
        using omni_type     = omni::stencil_t<
                grid::face_centered,
                omni::elem_t<omni::offset_t<-3, 0, 0>, info_type>,
                omni::elem_t<omni::offset_t<-1, 0, 0>, info_type>,
                omni::elem_t<omni::offset_t< 1, 0, 0>, info_type>,
                omni::elem_t<omni::offset_t< 3, 0, 0>, info_type>
            >;

        const flux_func_t flux_func;
        
        weno_t(const flux_func_t& flux_func_in)
        : flux_func{flux_func_in} {}
        
        _sp_hybrid output_type operator()(const auto& input) const
        {
            fluid_state::flux_t<float_t> output;
            auto fp0 = flux_func(input.cell(0_c));
            auto fp1 = flux_func(input.cell(1_c));
            auto fp2 = flux_func(input.cell(2_c));
            auto fp3 = flux_func(input.cell(3_c));
            const auto& f0_u = fp0[0];
            const auto& f0_d = fp0[1];
            const auto& f1_u = fp1[0];
            const auto& f1_d = fp1[1];
            const auto& f2_u = fp2[0];
            const auto& f2_d = fp2[1];
            const auto& f3_u = fp3[0];
            const auto& f3_d = fp3[1];
            for (std::size_t k = 0; k < output.size(); ++k)
            {
                //upwind
                const float_t r0  = -0.5*f0_u[k]+1.5*f1_u[k]; //sten0
                const float_t r1  =  0.5*f1_u[k]+0.5*f2_u[k]; //sten1
                
                //downwind
                const float_t r2  =  0.5*f1_d[k]+0.5*f2_d[k]; //sten2
                const float_t r3  =  1.5*f2_d[k]-0.5*f3_d[k]; //sten3
                
                const float_t be0 = (f0_u[k]-f1_u[k])*(f0_u[k]-f1_u[k]);
                const float_t be1 = (f1_u[k]-f2_u[k])*(f1_u[k]-f2_u[k]);
                
                const float_t be2 = (f1_d[k]-f2_d[k])*(f1_d[k]-f2_d[k]);
                const float_t be3 = (f2_d[k]-f3_d[k])*(f2_d[k]-f3_d[k]);
                
                const float_t eps = 1e-16;
                const float_t a0  = (1.0/3.0)/((be0+eps)*(be0+eps));
                const float_t a1  = (2.0/3.0)/((be1+eps)*(be1+eps));
                const float_t a2  = (2.0/3.0)/((be2+eps)*(be2+eps));
                const float_t a3  = (1.0/3.0)/((be3+eps)*(be3+eps));
                // const float_t a0  = (1.0/3.0)/((be0+eps));
                // const float_t a1  = (2.0/3.0)/((be1+eps));
                // const float_t a2  = (2.0/3.0)/((be2+eps));
                // const float_t a3  = (1.0/3.0)/((be3+eps));
                
                if constexpr (use_smooth == disable_smooth)
                {
                    //use this for MMS!!
                    output[k] = (1.0/3.0)*r0 + (2.0/3.0)*r1 + (2.0/3.0)*r2 + (1.0/3.0)*r3;
                }
                else
                {
                    const float_t w0 = a0/(a0+a1);
                    const float_t w1 = a1/(a0+a1);
                    const float_t w2 = a2/(a2+a3);
                    const float_t w3 = a3/(a2+a3);
                    output[k] = w0*r0 + w1*r1 + w2*r2 + w3*r3;
                }
            }
            
            return output;
        }
    };
}