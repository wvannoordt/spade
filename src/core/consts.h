#pragma once

#include "core/typedef.h"

namespace spade::consts
{
    constexpr static real_t pi = 3.1415926535;
    constexpr static real_t Rgas_uni = 8314.462;
	constexpr static real_t Na_kmol = 6.022045E26; // Avogadro's number in kilo moles
	constexpr static real_t Na = 6.022045E23; // Avogadro's number
	constexpr static real_t sigma = 5.6703744192E-8; // Stefan-Boltzmann constant for radiation
	constexpr static real_t kSI = 1.38064852E-23; // Boltzmann constant in SI units
	constexpr static real_t kCGS = 1.38064852E-16; // Boltzmann constant in CGS units
	constexpr static real_t eCGS = 4.8032E-10; // Electron charge in CGS system
}
