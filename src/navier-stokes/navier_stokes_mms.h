#pragma once

#include "navier-stokes/fluid_state.h"

namespace cvdf::navier_stokes_mms
{
    template <fluid_state::state_independent_gas gas_t, viscous_laws::state_independent_viscosity visc_t> struct cns_perfectgas_mms_t
    {
        typedef typename gas_t::value_type dtype; 
        gas_t gas;
        visc_t visc;
        const dtype pi = 3.1415926535;
        const dtype two_pi = 2.0*pi;
        
        cns_perfectgas_mms_t(const gas_t& gas_in, const visc_t& vlaw_in)
        {
            gas  = gas_in;
            visc = vlaw_in;
        }
        
        //P and derivatives
        dtype  P(const dtype& x, const dtype& y, const dtype& z)    const {return 5.0+2.0*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)+sin(4.0*two_pi*z);}
        dtype dP_dx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*3.0*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y);}
        dtype dP_dy(const dtype& x, const dtype& y, const dtype& z) const {return 2.0*2.0*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y);}
        dtype dP_dz(const dtype& x, const dtype& y, const dtype& z) const {return 4.0*two_pi*cos(4.0*two_pi*z);}

        //T and derivatives
        dtype  T(const dtype& x, const dtype& y, const dtype& z)    const {return 10.0+2.0*cos(2.0*two_pi*x)*sin(3.0*two_pi*y)+sin(4.0*two_pi*z);}
        dtype dT_dx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*sin(2.0*two_pi*x)*sin(3.0*two_pi*y);}
        dtype dT_dy(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*3.0*two_pi*cos(2.0*two_pi*x)*cos(3.0*two_pi*y);}
        dtype dT_dz(const dtype& x, const dtype& y, const dtype& z) const {return  4.0*two_pi*cos(4.0*two_pi*z);}
        dtype d2T_dx2(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*2.0*two_pi*two_pi*cos(2.0*two_pi*x)*sin(3.0*two_pi*y);}
        dtype d2T_dy2(const dtype& x, const dtype& y, const dtype& z) const {return  -2.0*3.0*3.0*two_pi*two_pi*cos(2.0*two_pi*x)*sin(3.0*two_pi*y);}
        dtype d2T_dz2(const dtype& x, const dtype& y, const dtype& z) const {return  -4.0*4.0*two_pi*two_pi*sin(4.0*two_pi*z);}
        
        //U and derivatives
        dtype  U(const dtype& x, const dtype& y, const dtype& z)    const {return sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype dU_dx(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype dU_dy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype dU_dz(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        
        dtype d2U_dxx(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dxy(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dxz(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        
        dtype d2U_dyx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dyy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dyz(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        
        dtype d2U_dzx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        dtype d2U_dzy(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        dtype d2U_dzz(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        
        //V and derivatives
        dtype  V(const dtype& x, const dtype& y, const dtype& z)    const {return cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype dV_dx(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype dV_dy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype dV_dz(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        
        dtype d2V_dxx(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dxy(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dxz(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        
        dtype d2V_dyx(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dyy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dyz(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        
        dtype d2V_dzx(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        dtype d2V_dzy(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        dtype d2V_dzz(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        
        //W and derivatives
        dtype  W(const dtype& x, const dtype& y, const dtype& z)    const {return sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype dW_dx(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype dW_dy(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype dW_dz(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        
        dtype d2W_dxx(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dxy(const dtype& x, const dtype& y, const dtype& z) const {return   3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dxz(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*4.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        
        dtype d2W_dyx(const dtype& x, const dtype& y, const dtype& z) const {return   2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dyy(const dtype& x, const dtype& y, const dtype& z) const {return  -2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dyz(const dtype& x, const dtype& y, const dtype& z) const {return  -2.0*4.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        
        dtype d2W_dzx(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        dtype d2W_dzy(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        dtype d2W_dzz(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*4.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        
        dtype  Rho(const dtype& x, const dtype& y, const dtype& z) const {return P(x,y,z)/(gas.get_R()*T(x,y,z));}
        dtype dRho_dx(const dtype& x, const dtype& y, const dtype& z) const
        {
            return (T(x,y,z)*dP_dx(x,y,z)-P(x,y,z)*dT_dx(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.get_R());
        }
        
        dtype dRho_dy(const dtype& x, const dtype& y, const dtype& z) const
        {
            return (T(x,y,z)*dP_dy(x,y,z)-P(x,y,z)*dT_dy(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.get_R());
        }
        
        dtype dRho_dz(const dtype& x, const dtype& y, const dtype& z) const
        {
            return (T(x,y,z)*dP_dz(x,y,z)-P(x,y,z)*dT_dz(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.get_R());
        }
        
        dtype  H(const dtype& x, const dtype& y, const dtype& z) const
        {
            return 0.5*(U(x,y,z)*U(x,y,z) + V(x,y,z)*V(x,y,z) + W(x,y,z)*W(x,y,z)) + gas.get_gamma()*P(x,y,z)/((gas.get_gamma()-1.0)*Rho(x,y,z));
        }
        dtype dH_dx(const dtype& x, const dtype& y, const dtype& z) const
        {
            return U(x,y,z)*dU_dx(x,y,z)+V(x,y,z)*dV_dx(x,y,z)+W(x,y,z)*dW_dx(x,y,z)
                + gas.get_gamma()*(Rho(x,y,z)*dP_dx(x,y,z)-dRho_dx(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.get_gamma()-1.0));
        }
        dtype dH_dy(const dtype& x, const dtype& y, const dtype& z) const
        {
            return U(x,y,z)*dU_dy(x,y,z)+V(x,y,z)*dV_dy(x,y,z)+W(x,y,z)*dW_dy(x,y,z)
                + gas.get_gamma()*(Rho(x,y,z)*dP_dy(x,y,z)-dRho_dy(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.get_gamma()-1.0));
        }
        dtype dH_dz(const dtype& x, const dtype& y, const dtype& z) const
        {
            return U(x,y,z)*dU_dz(x,y,z)+V(x,y,z)*dV_dz(x,y,z)+W(x,y,z)*dW_dz(x,y,z)
                + gas.get_gamma()*(Rho(x,y,z)*dP_dz(x,y,z)-dRho_dz(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.get_gamma()-1.0));
        }
        
        ctrs::array<dtype,5> conv_rhs(const ctrs::array<dtype,3>& x) const
        {
            ctrs::array<dtype,5> rhsAr;
            rhsAr[0] = 
                -U(x[0],x[1],x[2])*dRho_dx(x[0],x[1],x[2])-dU_dx(x[0],x[1],x[2])*Rho(x[0],x[1],x[2])
                -V(x[0],x[1],x[2])*dRho_dy(x[0],x[1],x[2])-dV_dy(x[0],x[1],x[2])*Rho(x[0],x[1],x[2])
                -W(x[0],x[1],x[2])*dRho_dz(x[0],x[1],x[2])-dW_dz(x[0],x[1],x[2])*Rho(x[0],x[1],x[2]);
            rhsAr[1] = 
                -dRho_dx(x[0],x[1],x[2])*U    (x[0],x[1],x[2])*H    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dU_dx(x[0],x[1],x[2])*H    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])*dH_dx(x[0],x[1],x[2])
                -dRho_dy(x[0],x[1],x[2])*V    (x[0],x[1],x[2])*H    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dV_dy(x[0],x[1],x[2])*H    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])*dH_dy(x[0],x[1],x[2])
                -dRho_dz(x[0],x[1],x[2])*W    (x[0],x[1],x[2])*H    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dW_dz(x[0],x[1],x[2])*H    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])*dH_dz(x[0],x[1],x[2]);
            rhsAr[2] = 
                -dRho_dx(x[0],x[1],x[2])*U    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dU_dx(x[0],x[1],x[2])*U    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])*dU_dx(x[0],x[1],x[2])
                -dRho_dy(x[0],x[1],x[2])*V    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dV_dy(x[0],x[1],x[2])*U    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])*dU_dy(x[0],x[1],x[2])
                -dRho_dz(x[0],x[1],x[2])*W    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dW_dz(x[0],x[1],x[2])*U    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])*dU_dz(x[0],x[1],x[2])
                -dP_dx  (x[0],x[1],x[2]);
            rhsAr[3] = 
                -dRho_dx(x[0],x[1],x[2])*U    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dU_dx(x[0],x[1],x[2])*V    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])*dV_dx(x[0],x[1],x[2])
                -dRho_dy(x[0],x[1],x[2])*V    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dV_dy(x[0],x[1],x[2])*V    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])*dV_dy(x[0],x[1],x[2])
                -dRho_dz(x[0],x[1],x[2])*W    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dW_dz(x[0],x[1],x[2])*V    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])*dV_dz(x[0],x[1],x[2])
                -dP_dy  (x[0],x[1],x[2]);
            rhsAr[4] = 
                -dRho_dx(x[0],x[1],x[2])*U    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dU_dx(x[0],x[1],x[2])*W    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*U    (x[0],x[1],x[2])*dW_dx(x[0],x[1],x[2])
                -dRho_dy(x[0],x[1],x[2])*V    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dV_dy(x[0],x[1],x[2])*W    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*V    (x[0],x[1],x[2])*dW_dy(x[0],x[1],x[2])
                -dRho_dz(x[0],x[1],x[2])*W    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*dW_dz(x[0],x[1],x[2])*W    (x[0],x[1],x[2])
                -Rho    (x[0],x[1],x[2])*W    (x[0],x[1],x[2])*dW_dz(x[0],x[1],x[2])
                -dP_dz  (x[0],x[1],x[2]);
            return rhsAr;
        }
        
        dtype Tau_xx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dU_dx(x,y,z)+dU_dx(x,y,z))+visc.get_beta()*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
        dtype Tau_xy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dU_dy(x,y,z)+dV_dx(x,y,z)); }
        dtype Tau_xz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dU_dz(x,y,z)+dW_dx(x,y,z)); }
        dtype Tau_yx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dV_dx(x,y,z)+dU_dy(x,y,z)); }
        dtype Tau_yy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dV_dy(x,y,z)+dV_dy(x,y,z))+visc.get_beta()*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
        dtype Tau_yz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dV_dz(x,y,z)+dW_dy(x,y,z)); }
        dtype Tau_zx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dW_dx(x,y,z)+dU_dz(x,y,z)); }
        dtype Tau_zy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dW_dy(x,y,z)+dV_dz(x,y,z)); }
        dtype Tau_zz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(dW_dz(x,y,z)+dW_dz(x,y,z))+visc.get_beta()*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
        
        dtype dTau_xx_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dxx(x,y,z)+d2U_dxx(x,y,z))+visc.get_beta()*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
        dtype dTau_xx_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dxy(x,y,z)+d2U_dxy(x,y,z))+visc.get_beta()*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
        dtype dTau_xx_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dxz(x,y,z)+d2U_dxz(x,y,z))+visc.get_beta()*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
        
        dtype dTau_xy_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dyx(x,y,z)+d2V_dxx(x,y,z)); }
        dtype dTau_xy_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dyy(x,y,z)+d2V_dxy(x,y,z)); }
        dtype dTau_xy_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dyz(x,y,z)+d2V_dxz(x,y,z)); }
        
        dtype dTau_xz_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dzx(x,y,z)+d2W_dxx(x,y,z)); }
        dtype dTau_xz_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dzy(x,y,z)+d2W_dxy(x,y,z)); }
        dtype dTau_xz_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2U_dzz(x,y,z)+d2W_dxz(x,y,z)); }
        
        dtype dTau_yx_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dxx(x,y,z)+d2U_dyx(x,y,z)); }
        dtype dTau_yx_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dxy(x,y,z)+d2U_dyy(x,y,z)); }
        dtype dTau_yx_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dxz(x,y,z)+d2U_dyz(x,y,z)); }
        
        dtype dTau_yy_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dyx(x,y,z)+d2V_dyx(x,y,z))+visc.get_beta()*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
        dtype dTau_yy_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dyy(x,y,z)+d2V_dyy(x,y,z))+visc.get_beta()*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
        dtype dTau_yy_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dyz(x,y,z)+d2V_dyz(x,y,z))+visc.get_beta()*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
        
        dtype dTau_yz_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dzx(x,y,z)+d2W_dyx(x,y,z)); }
        dtype dTau_yz_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dzy(x,y,z)+d2W_dyy(x,y,z)); }
        dtype dTau_yz_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2V_dzz(x,y,z)+d2W_dyz(x,y,z)); }
        
        dtype dTau_zx_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dxx(x,y,z)+d2U_dzx(x,y,z)); }
        dtype dTau_zx_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dxy(x,y,z)+d2U_dzy(x,y,z)); }
        dtype dTau_zx_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dxz(x,y,z)+d2U_dzz(x,y,z)); }
        
        dtype dTau_zy_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dyx(x,y,z)+d2V_dzx(x,y,z)); }
        dtype dTau_zy_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dyy(x,y,z)+d2V_dzy(x,y,z)); }
        dtype dTau_zy_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dyz(x,y,z)+d2V_dzz(x,y,z)); }
        
        dtype dTau_zz_dx(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dzx(x,y,z)+d2W_dzx(x,y,z))+visc.get_beta()*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
        dtype dTau_zz_dy(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dzy(x,y,z)+d2W_dzy(x,y,z))+visc.get_beta()*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
        dtype dTau_zz_dz(const dtype& x, const dtype& y, const dtype& z) const { return visc.get_visc()*(d2W_dzz(x,y,z)+d2W_dzz(x,y,z))+visc.get_beta()*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
        
        dtype dQ_x_dx(const dtype& x, const dtype& y, const dtype& z) const { return (visc.get_conductivity())*d2T_dx2(x,y,z); }
        dtype dQ_y_dy(const dtype& x, const dtype& y, const dtype& z) const { return (visc.get_conductivity())*d2T_dy2(x,y,z); }
        dtype dQ_z_dz(const dtype& x, const dtype& y, const dtype& z) const { return (visc.get_conductivity())*d2T_dz2(x,y,z); }
        
        ctrs::array<dtype,5> visc_rhs(const ctrs::array<dtype,3>& x) const
        {
            ctrs::array<dtype,5> rhsAr;
            rhsAr[0] = 0.0;
            rhsAr[1] = 
                U(x[0],x[1],x[2])*dTau_xx_dx(x[0],x[1],x[2])+dU_dx(x[0],x[1],x[2])*Tau_xx(x[0],x[1],x[2])+
                V(x[0],x[1],x[2])*dTau_xy_dx(x[0],x[1],x[2])+dV_dx(x[0],x[1],x[2])*Tau_xy(x[0],x[1],x[2])+
                W(x[0],x[1],x[2])*dTau_xz_dx(x[0],x[1],x[2])+dW_dx(x[0],x[1],x[2])*Tau_xz(x[0],x[1],x[2])+
                U(x[0],x[1],x[2])*dTau_yx_dy(x[0],x[1],x[2])+dU_dy(x[0],x[1],x[2])*Tau_yx(x[0],x[1],x[2])+
                V(x[0],x[1],x[2])*dTau_yy_dy(x[0],x[1],x[2])+dV_dy(x[0],x[1],x[2])*Tau_yy(x[0],x[1],x[2])+
                W(x[0],x[1],x[2])*dTau_yz_dy(x[0],x[1],x[2])+dW_dy(x[0],x[1],x[2])*Tau_yz(x[0],x[1],x[2])+
                U(x[0],x[1],x[2])*dTau_zx_dz(x[0],x[1],x[2])+dU_dz(x[0],x[1],x[2])*Tau_zx(x[0],x[1],x[2])+
                V(x[0],x[1],x[2])*dTau_zy_dz(x[0],x[1],x[2])+dV_dz(x[0],x[1],x[2])*Tau_zy(x[0],x[1],x[2])+
                W(x[0],x[1],x[2])*dTau_zz_dz(x[0],x[1],x[2])+dW_dz(x[0],x[1],x[2])*Tau_zz(x[0],x[1],x[2])-
                dQ_x_dx(x[0],x[1],x[2])-dQ_y_dy(x[0],x[1],x[2])-dQ_z_dz(x[0],x[1],x[2]);
            rhsAr[2] = dTau_xx_dx(x[0],x[1],x[2])+dTau_xy_dy(x[0],x[1],x[2])+dTau_xz_dz(x[0],x[1],x[2]);
            rhsAr[3] = dTau_yx_dx(x[0],x[1],x[2])+dTau_yy_dy(x[0],x[1],x[2])+dTau_yz_dz(x[0],x[1],x[2]);
            rhsAr[4] = dTau_zx_dx(x[0],x[1],x[2])+dTau_zy_dy(x[0],x[1],x[2])+dTau_zz_dz(x[0],x[1],x[2]);
            return rhsAr;
        }
        
        fluid_state::prim_t<dtype> test_fcn(const ctrs::array<dtype,3>& x) const
        {
            fluid_state::prim_t<dtype> primsAr;
            primsAr[0] = P(x[0], x[1], x[2]);
            primsAr[1] = T(x[0], x[1], x[2]);
            primsAr[2] = U(x[0], x[1], x[2]);
            primsAr[3] = V(x[0], x[1], x[2]);
            primsAr[4] = W(x[0], x[1], x[2]);
            return primsAr;
        }
    };
}