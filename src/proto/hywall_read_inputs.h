#pragma once

#include "PTL.h"
#include "HyWall.h"

namespace spade::proto
{
    bool streq(const std::string& s1, const std::string& s2) {return s1.compare(s2)==0;}
    static void hywall_read_inputs(PTL::PropertySection& input_section, HyWall::UserSettings& settings)
    {
        settings.enableWallModel = true;
        settings.readRestart = false;
        input_section["rayDim"].MapTo(&settings.rayDim)                      = new PTL::PTLInteger(30, "number of ray points");
        input_section["solveSkip"].MapTo(&settings.solveSkip)                      = new PTL::PTLInteger(1, "solve interval");
        settings.asyncSolve = false;
        input_section["verboseLevel"].MapTo(&settings.verboseLevel)          = new PTL::PTLInteger(1, "debug output level");
        input_section["maxIterations"].MapTo(&settings.maxIterations)        = new PTL::PTLInteger(100, "Max. iterations");
        input_section["wallSpacing"].MapTo(&settings.wallSpacing)            = new PTL::PTLDouble(1e-6, "Max. itenumber of ray pointsrations");
        input_section["wallTemperature"].MapTo(&settings.wallTemperature)    = new PTL::PTLDouble(100, "Wall Temperature");
        input_section["adiabaticWall"].MapTo(&settings.adiabaticWall)        = new PTL::PTLBoolean(true, "Adiabatic wall");
        input_section["variablePrandtlT"].MapTo(&settings.variablePrandtlT)  = new PTL::PTLBoolean(false, "Variable turbulent prandtl number");
        input_section["fluidCp"].MapTo(&settings.fluidCp)                    = new PTL::PTLDouble(1005.0, "Specific heat");
        input_section["turbPradntl"].MapTo(&settings.turbPradntl)            = new PTL::PTLDouble(0.72, "Turbulent Prandtl");
        input_section["fluidPrandtl"].MapTo(&settings.fluidPrandtl)          = new PTL::PTLDouble(0.9, "Laminar Prandtl");
        input_section["vanDriestAPlus"].MapTo(&settings.vanDriestAPlus)      = new PTL::PTLDouble(17.0, "van Driest Constant");
        input_section["gasConstant"].MapTo(&settings.gasConstant)            = new PTL::PTLDouble(287.0, "Gas constant");
        input_section["enableTransitionSensor"].MapTo(&settings.enableTransitionSensor) = new PTL::PTLBoolean(false, "Enable Transition Sensor");
        
        std::string mom_eq_str, trb_eq_str, eng_eq_str;
        input_section["momentumEquationType"].MapTo(&mom_eq_str)          = new PTL::PTLString("ODE", "Momentum equation type");
        input_section["turbulenceEquationType"].MapTo(&trb_eq_str)        = new PTL::PTLString("vanDriest", "Turbulence equation type");
        input_section["energyEquationType"].MapTo(&eng_eq_str)            = new PTL::PTLString("ODE", "Energy equation type");
        
        input_section["momentumUnderRelaxationODE"].MapTo(&settings.momentumUnderRelaxationODE)     = new PTL::PTLDouble(0.8, "Relaxation factor for momentum ODE");
        input_section["turbulenceUnderRelaxationODE"].MapTo(&settings.turbulenceUnderRelaxationODE) = new PTL::PTLDouble(0.6, "Relaxation factor for turbulence ODE");
        input_section["energyUnderRelaxationODE"].MapTo(&settings.energyUnderRelaxationODE)         = new PTL::PTLDouble(0.7, "Relaxation factor for energy ODE");
        input_section["includeMomentumRhs"].MapTo(&settings.includeMomentumRhs)                     = new PTL::PTLBoolean(false, "Include the parameterized convection term");
        input_section["isCompressible"].MapTo(&settings.isCompressible)                             = new PTL::PTLBoolean(false, "Use variable density");
        input_section["suthViscRef"].MapTo(&settings.suthViscRef)                                   = new PTL::PTLDouble(1.45151376745308e-06, "Reference viscosity for viscosity power law");
        input_section["suthTRef"].MapTo(&settings.suthTRef)                                         = new PTL::PTLDouble(110.4, "Reference temperature for viscosity power law");
        
        std::string visc_law_str;
        input_section["viscousLaw"].MapTo(&visc_law_str)                                            = new PTL::PTLString("constant", "Viscous law");
        
        std::string yscale_str;
        input_section["yScale"].MapTo(&yscale_str) = new PTL::PTLString("trettelLarsson", "y-coordinate scaling");
        
        bool failed = false;
        if (!failed) input_section.StrictParse();
             if (streq(mom_eq_str, "allmaras")) {settings.momentumEquationType = HyCore::momentum::allmaras;}
        else if (streq(mom_eq_str, "ODE"))      {settings.momentumEquationType = HyCore::momentum::ODE;}
        else if (!failed) {std::cout << "Invalid value for momentum equation" << std::endl; abort();}
        
             if (streq(trb_eq_str, "linear"))    {settings.turbulenceEquationType = HyCore::turbulence::linear;}
        else if (streq(trb_eq_str, "ODE"))       {settings.turbulenceEquationType = HyCore::turbulence::ODE;}
        else if (streq(trb_eq_str, "vanDriest")) {settings.turbulenceEquationType = HyCore::turbulence::vanDriest;}
        else if (!failed) {std::cout << "Invalid value for turbulence equation" << std::endl; abort();}
        
             if (streq(eng_eq_str, "croccoBusemann")) {settings.energyEquationType = HyCore::energy::croccoBusemann;}
        else if (streq(eng_eq_str, "ODE"))            {settings.energyEquationType = HyCore::energy::ODE; }
        else if (streq(eng_eq_str, "linear"))         {settings.energyEquationType = HyCore::energy::linear;}
        else if (!failed) {std::cout << "Invalid value for energy equation" << std::endl; abort();}
        
            if  (streq(yscale_str, "trettelLarsson")) {settings.yscaleType = HyCore::yscale::trettelLarsson;}
        else if (streq(yscale_str, "yPlus"))          {settings.yscaleType = HyCore::yscale::yPlus;}
        else if (streq(yscale_str, "mixed"))          {settings.yscaleType = HyCore::yscale::mixed;}
        else if (!failed) {std::cout << "Invalid value for y-coordinate transform" << std::endl; abort();}
        
             if (streq(visc_law_str, "constant"))   {settings.viscousLaw = HyCore::visclaw::constant;}
        else if (streq(visc_law_str, "sutherland")) {settings.viscousLaw = HyCore::visclaw::sutherland;}
        else if (streq(visc_law_str, "PowerLaw"))   {settings.viscousLaw = HyCore::visclaw::PowerLaw;}
        else if (!failed) {std::cout << "Invalid value for energy equation" << std::endl; abort();}
    }
}