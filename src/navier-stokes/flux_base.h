#pragma once

namespace cvdf::flux_base
{
    template <typename d> struct flux_input_t
    {
        
    }
    template <class derived_t, class output_t, class input_t> struct flux_base_t
    {
        virtual output_t calc_flux(const input_t& flux_input) const = 0;
    };
}