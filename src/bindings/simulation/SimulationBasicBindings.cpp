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
namespace simulation {

void init_basic_bindings(py::module& m) {
    // Check if LogLevel is already defined (it might be defined elsewhere)
    if (!py::hasattr(m, "LogLevel")) {
        // System logging bindings (using new SystemLogger)
        py::enum_<system::common::LogLevel>(m, "LogLevel")
            .value("DEBUG", system::common::LogLevel::DEBUG)
            .value("INFO", system::common::LogLevel::INFO)
            .value("WARNING", system::common::LogLevel::WARNING)
            .value("ERROR", system::common::LogLevel::ERROR);
    }
    
    // Platform logging (keep for compatibility, with different name)
    if (!py::hasattr(m, "PlatformLogLevel")) {
        py::enum_<platform::LogLevel>(m, "PlatformLogLevel")
            .value("DEBUG", platform::LogLevel::DEBUG)
            .value("INFO", platform::LogLevel::INFO)
            .value("WARNING", platform::LogLevel::WARNING)
            .value("ERROR", platform::LogLevel::ERROR);
    }
    
    // System logger functions
    m.def("set_verbose", &system::common::SystemLogger::setVerbose,
         "Set verbose mode for system logging");
    m.def("set_log_level", 
         [](system::common::LogLevel level) { system::common::SystemLogger::setLogLevel(level); },
         "Set the minimum log level for system");
    m.def("set_debug_mode", &system::common::SystemLogger::setDebugMode,
         "Enable or disable debug mode for system logging");
    
    // Platform logging (keep for compatibility)
    m.def("set_platform_verbose", &platform::set_verbose, 
         "Set verbose mode for platform logging");
    m.def("set_platform_log_level", &platform::set_log_level,
         "Set the minimum log level for platform");
    m.def("set_platform_debug_mode", &platform::set_debug_mode,
         "Enable or disable debug mode for platform logging");

    // Basic energy calculation functions (using new EnergyAPI)
    m.def("computeMovementEnergy", &platform::cpu::energy::computeMovementEnergy,
          "Calculate nonbonded energies for movement residues only");
          
    m.def("computeMovementEnergyCutoff", &platform::cpu::energy::computeMovementEnergyCutoff,
          "Calculate nonbonded energies for movement residues only with distance cutoff");
          
    // Keep old functions for backward compatibility (they return None)
    m.def("_computeSystemEnergyVoid", &platform::cpu::energy::computeSystemEnergy,
          "Calculate nonbonded energies for the full system (internal, updates state)");
          
    m.def("_computeSystemEnergyCutoffVoid", &platform::cpu::energy::computeSystemEnergyCutoff,
          "Calculate nonbonded energies for the full system with distance cutoff (internal, updates state)");
    
    // New functions that return total energy  
    m.def("computeSystemEnergy", &platform::cpu::energy::computeSystemEnergyTotal,
          "Calculate and return total system energy");
          
    m.def("computeSystemEnergyCutoff", &platform::cpu::energy::computeSystemEnergyCutoffTotal,
          "Calculate and return total system energy with distance cutoff");
          
    m.def("computeSystemEnergyPBC", &platform::cpu::energy::computeSystemEnergyPBC,
          "Calculate nonbonded energies for the full system with periodic boundary conditions");
          
    m.def("computeSystemEnergyPBCCutoff", &platform::cpu::energy::computeSystemEnergyPBCCutoff,
          "Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff");
    
    m.def("computeSystemVdwEnergyCutoff", &platform::cpu::energy::computeSystemVdwEnergyCutoff,
          "Calculate VDW energies for the full system with distance cutoff");
          
    m.def("setEnergyDebugOutput", &platform::cpu::energy::setEnergyDebugOutput,
          "Enable or disable debug output for energy calculations");
    
    m.def("getTotalEnergyComponents", &platform::cpu::energy::getTotalEnergyComponents,
          "Get the total electrostatic and van der Waals energy components as a tuple (elec, vdw)");
    
    // CHARMM switching function related bindings have been removed
    // Please use the set_switching_function and calculate_switching_function methods in the MonteCarloSystem class
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc