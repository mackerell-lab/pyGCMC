// src/bindings/simulation/SimulationBasicBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../simulation/simulation.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_basic_bindings(py::module& m) {
    // Platform logging bindings
    py::enum_<platform::LogLevel>(m, "PlatformLogLevel")
        .value("DEBUG", platform::LogLevel::DEBUG)
        .value("INFO", platform::LogLevel::INFO)
        .value("WARNING", platform::LogLevel::WARNING)
        .value("ERROR", platform::LogLevel::ERROR);
    
    m.def("set_platform_verbose", &platform::set_verbose, 
         "Set verbose mode for platform logging");
    m.def("set_platform_log_level", &platform::set_log_level,
         "Set the minimum log level for platform");
    m.def("set_platform_debug_mode", &::pygcmc::simulation::set_debug_mode,
         "Enable or disable debug mode for platform logging");

    // Basic energy calculation functions
    m.def("computeMovementEnergy", &::pygcmc::simulation::Simulation::computeMovementEnergy,
          "Calculate nonbonded energies for movement residues only");
          
    m.def("computeMovementEnergyCutoff", &::pygcmc::simulation::Simulation::computeMovementEnergyCutoff,
          "Calculate nonbonded energies for movement residues only with distance cutoff");
          
    m.def("computeSystemEnergy", &::pygcmc::simulation::Simulation::computeSystemEnergy,
          "Calculate nonbonded energies for the full system");
          
    m.def("computeSystemEnergyCutoff", &::pygcmc::simulation::Simulation::computeSystemEnergyCutoff,
          "Calculate nonbonded energies for the full system with distance cutoff");
          
    m.def("computeSystemEnergyPBC", &::pygcmc::simulation::Simulation::computeSystemEnergyPBC,
          "Calculate nonbonded energies for the full system with periodic boundary conditions");
          
    m.def("computeSystemEnergyPBCCutoff", &::pygcmc::simulation::Simulation::computeSystemEnergyPBCCutoff,
          "Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff");
    
    m.def("computeSystemVdwEnergyCutoff", &::pygcmc::simulation::Simulation::computeSystemVdwEnergyCutoff,
          "Calculate VDW energies for the full system with distance cutoff");
          
    m.def("setEnergyDebugOutput", &::pygcmc::simulation::Simulation::setEnergyDebugOutput,
          "Enable or disable debug output for energy calculations");
    
    // CHARMM switching function related bindings have been removed
    // Please use the set_switching_function and calculate_switching_function methods in the MonteCarloSystem class
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc