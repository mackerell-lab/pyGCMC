// GCMCBindings.cpp - Simplified Python bindings for GCMC simulation

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// Use new EnergyAPI instead of simulation.hpp
#include "../../platform/cpu/energy/EnergyAPI.hpp"
#include "../../platform/cpu/energy/common/EnergyInterface.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_gcmc_bindings(py::module& m) {
    // Energy method enum
    py::enum_<platform::cpu::EnergyMethod>(m, "EnergyMethod")
        .value("DIRECT", platform::cpu::EnergyMethod::DIRECT)
        .value("EWALD", platform::cpu::EnergyMethod::EWALD)
        .value("PME", platform::cpu::EnergyMethod::PME);
    
    // For now, we skip the GCMCSimulation class bindings since
    // that functionality has been moved to MovementAPI
    // The actual GCMC functionality is exposed through movement module
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc