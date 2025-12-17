// src/bindings/simulation/SimulationBasicBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// Use new API headers instead of simulation.hpp
#include "../../platform/cpu/energy/EnergyAPI.hpp"
#include "../../system/common/SystemLogger.hpp"
#include "../../platform/platform.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace platform {

void init_basic_bindings(py::module& m) {
    // Check if LogLevel is already defined (it might be defined elsewhere)
    if (!py::hasattr(m, "LogLevel")) {
        // System logging bindings (using new SystemLogger)
        py::enum_<::pygcmc::system::common::LogLevel>(m, "LogLevel")
            .value("DEBUG", ::pygcmc::system::common::LogLevel::DEBUG)
            .value("INFO", ::pygcmc::system::common::LogLevel::INFO)
            .value("WARNING", ::pygcmc::system::common::LogLevel::WARNING)
            .value("ERROR", ::pygcmc::system::common::LogLevel::ERROR);
    }
    
    // Platform logging (keep for compatibility, with different name)
    if (!py::hasattr(m, "PlatformLogLevel")) {
        py::enum_<::pygcmc::platform::LogLevel>(m, "PlatformLogLevel")
            .value("DEBUG", ::pygcmc::platform::LogLevel::DEBUG)
            .value("INFO", ::pygcmc::platform::LogLevel::INFO)
            .value("WARNING", ::pygcmc::platform::LogLevel::WARNING)
            .value("ERROR", ::pygcmc::platform::LogLevel::ERROR);
    }
    
    // System logger functions
    m.def("set_verbose", &::pygcmc::system::common::SystemLogger::setVerbose,
         "Set verbose mode for system logging");
    m.def("set_log_level", 
         [](::pygcmc::system::common::LogLevel level) { ::pygcmc::system::common::SystemLogger::setLogLevel(level); },
         "Set the minimum log level for system");
    m.def("set_debug_mode", &::pygcmc::system::common::SystemLogger::setDebugMode,
         "Enable or disable debug mode for system logging");
    
    // Platform logging (keep for compatibility)
    m.def("set_platform_verbose", &::pygcmc::platform::set_verbose, 
         "Set verbose mode for platform logging");
    m.def("set_platform_log_level", &::pygcmc::platform::set_log_level,
         "Set the minimum log level for platform");
    m.def("set_platform_debug_mode", &::pygcmc::platform::set_debug_mode,
         "Enable or disable debug mode for platform logging");

    // Basic energy calculation functions (using new EnergyAPI)
    m.def("computeMovementEnergy", &::pygcmc::platform::cpu::energy::computeMovementEnergy,
          "Calculate nonbonded energies for movement residues only");
          
    m.def("computeMovementEnergyCutoff", &::pygcmc::platform::cpu::energy::computeMovementEnergyCutoff,
          "Calculate nonbonded energies for movement residues only with distance cutoff");
          
    // Keep old functions for backward compatibility (they return None)
    m.def("_computeSystemEnergyVoid", &::pygcmc::platform::cpu::energy::computeSystemEnergy,
          "Calculate nonbonded energies for the full system (internal, updates state)");
          
    m.def("_computeSystemEnergyCutoffVoid", &::pygcmc::platform::cpu::energy::computeSystemEnergyCutoff,
          "Calculate nonbonded energies for the full system with distance cutoff (internal, updates state)");
    
    // New functions that return total energy  
    m.def("computeSystemEnergy", &::pygcmc::platform::cpu::energy::computeSystemEnergyTotal,
          "Calculate and return total system energy");
          
    m.def("computeSystemEnergyCutoff", &::pygcmc::platform::cpu::energy::computeSystemEnergyCutoffTotal,
          "Calculate and return total system energy with distance cutoff");
          
    m.def("computeSystemEnergyPBC", &::pygcmc::platform::cpu::energy::computeSystemEnergyPBC,
          "Calculate nonbonded energies for the full system with periodic boundary conditions");
          
    m.def("computeSystemEnergyPBCCutoff", &::pygcmc::platform::cpu::energy::computeSystemEnergyPBCCutoff,
          "Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff");
    
    m.def("computeSystemVdwEnergyCutoff", &::pygcmc::platform::cpu::energy::computeSystemVdwEnergyCutoff,
          "Calculate VDW energies for the full system with distance cutoff");
          
    m.def("setEnergyDebugOutput", &::pygcmc::platform::cpu::energy::setEnergyDebugOutput,
          "Enable or disable debug output for energy calculations");
    
    m.def("getTotalEnergyComponents", &::pygcmc::platform::cpu::energy::getTotalEnergyComponents,
          "Get the total electrostatic and van der Waals energy components as a tuple (elec, vdw)");

    // Total energy helpers (do not recompute energies; uses current state fields)
    m.def("getTotalEnergyUniquePairs", &::pygcmc::platform::cpu::energy::getTotalEnergyUniquePairs,
          py::arg("state"), py::arg("method"),
          "Get total energy using a unique-pairs convention for pair interactions");
    
    // CHARMM switching function related bindings have been removed
    // Please use the set_switching_function and calculate_switching_function methods in the MonteCarloSystem class
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc
