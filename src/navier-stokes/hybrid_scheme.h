#pragma once

namespace spade::convective
{
    template <typename scheme0_t, typename scheme1_t, typename blender_t>
    struct hybrid_scheme_t
    {
        const scheme0_t& scheme0;
        const scheme1_t& scheme1;
        const blender_t& blender;
        
        hybrid_scheme_t(const scheme0_t& scheme0_in, const scheme1_t& scheme1_in, const blender_t& blender_in)
        : scheme0{scheme0_in}, scheme1{scheme1_in}, blender{blender_in} {}
    };
}