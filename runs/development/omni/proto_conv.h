#pragma once

//ignore the contents of this file

namespace proto
{
    template <spade::fluid_state::state_dependent_gas gas_t> struct totani_lr
    {
        using real_type         = gas_t::value_type;
        using output_type       = spade::fluid_state::flux_t<real_type>;
        
        // using face_stencil_type = spade::omni::directed_stencil_t
        // <
        //     spade::grid::face_centered,
        //     spade::omni::mono_t
        // >;
        
        // using cell_stencil_type = spade::omni::directed_stencil_t
        // <
        //     spade::grid::cell_centered,
        //     spade::omni::line_t<2>
        // >;
        
        // using stencil_type = spade::omni::union_t<>
        
        totani_lr(const gas_t& gas_in) {gas = &gas_in;}
        
        
        // output_type calc_flux(const input_type& input) const
        // {
        //     output_type output;
        //     const auto& ql       = std::get<0>(input.cell_data.left.elements).data;
        //     const auto& qr       = std::get<0>(input.cell_data.right.elements).data;
        //     const auto& normal_l = std::get<1>(input.cell_data.left.elements).data;
        //     const auto& normal_r = std::get<1>(input.cell_data.right.elements).data;
        // 
        //     real_type un_l = normal_l[0]*ql.u()+normal_l[1]*ql.v()+normal_l[2]*ql.w();
        //     real_type un_r = normal_r[0]*qr.u()+normal_r[1]*qr.v()+normal_r[2]*qr.w();
        // 
        //     real_type rho_l = ql.p()/(gas->get_R(ql)*ql.T());
        //     real_type rho_r = qr.p()/(gas->get_R(qr)*qr.T());
        // 
        //     real_type e_l = ql.p()/(rho_l*(gas->get_gamma(ql)-real_type(1.0)));
        //     real_type e_r = qr.p()/(rho_r*(gas->get_gamma(qr)-1.0));
        // 
        //     real_type c = 0.25*(rho_l+rho_r)*(un_l+un_r);
        //     output.continuity() = c;
        //     output.energy()     = 0.5*c*(e_l + e_r + ql.u()*qr.u() + ql.v()*qr.v() + ql.w()*qr.w()) + 0.5*(un_l*qr.p() + un_r*ql.p());
        //     output.x_momentum() = 0.5*c*(ql.u()+qr.u()) + 0.5*(normal_l[0]*ql.p()+normal_r[0]*qr.p());
        //     output.y_momentum() = 0.5*c*(ql.v()+qr.v()) + 0.5*(normal_l[1]*ql.p()+normal_r[1]*qr.p());
        //     output.z_momentum() = 0.5*c*(ql.w()+qr.w()) + 0.5*(normal_l[2]*ql.p()+normal_r[2]*qr.p());
        //     return output;
        // }
        
        const gas_t* gas;
    };
}