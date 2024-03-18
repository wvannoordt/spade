#pragma once

namespace spade::convective
{
    template <typename scheme0_t, typename scheme1_t, typename blender_t>
    struct hybrid_scheme_t
    {
        const scheme0_t scheme0; // centered
        const scheme1_t scheme1; // dissipative
        const blender_t blender;
        
        using output_type     = typename scheme0_t::output_type;
        using blend_stencil_t = omni::prefab::face_mono_t<typename blender_t::info_type>;
        using omni_type       = omni::stencil_union<blend_stencil_t, typename scheme0_t::omni_type, typename scheme1_t::omni_type>;
        
        _sp_hybrid hybrid_scheme_t(const scheme0_t& scheme0_in, const scheme1_t& scheme1_in, const blender_t& blender_in)
        : scheme0{scheme0_in}, scheme1{scheme1_in}, blender{blender_in} {}
        
        _sp_hybrid output_type operator() (const auto& input_data) const
        {
            auto blend = blender.get_sensor(omni::interpret_stencil<blend_stencil_t>(input_data).face(0_c));
            auto f0 = scheme0(omni::interpret_stencil<typename scheme0_t::omni_type>(input_data));
            auto f1 = scheme1(omni::interpret_stencil<typename scheme1_t::omni_type>(input_data));
            auto one = decltype(blend)(1.0);
            f0 *= (one-blend);
            f1 *= blend;
            f0 += f1;
            return f0;
        }
    };

    template <typename scheme0_t, typename scheme1_t, typename blender_t>
    struct hybrid_df_scheme_t
    {
        const scheme0_t scheme0; // centered
        const scheme1_t scheme1; // dissipative
        const blender_t blender;
        
        using output_type     = typename scheme0_t::output_type;
        using blend_stencil_t = omni::prefab::face_mono_t<typename blender_t::info_type>;
        using omni_type       = omni::stencil_union<blend_stencil_t, typename scheme0_t::omni_type, typename scheme1_t::omni_type>;
        
        _sp_hybrid hybrid_df_scheme_t(const scheme0_t& scheme0_in, const scheme1_t& scheme1_in, const blender_t& blender_in)
        : scheme0{scheme0_in}, scheme1{scheme1_in}, blender{blender_in} {}
        
        _sp_hybrid output_type operator() (const auto& input_data) const
        {
            auto blend = blender.get_sensor(omni::interpret_stencil<blend_stencil_t>(input_data).face(0_c));
            auto f0 = scheme0(omni::interpret_stencil<typename scheme0_t::omni_type>(input_data));
            auto f1 = scheme1(omni::interpret_stencil<typename scheme1_t::omni_type>(input_data));
            auto one = decltype(blend)(1.0);
            f1 *= blend;
            f0 += f1;
            return f0;
        }
    };
}
