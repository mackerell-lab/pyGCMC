// GCMCBindings.cpp - Simplified Python bindings for GCMC simulation

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// Use new EnergyAPI instead of simulation.hpp
#include "../../platform/cpu/energy/EnergyAPI.hpp"
#include "../../platform/cpu/energy/common/EnergyInterface.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace platform {

void init_gcmc_bindings(py::module& m) {
    // Energy method enum
    py::enum_<::pygcmc::platform::cpu::EnergyMethod>(m, "EnergyMethod")
        .value("DIRECT", ::pygcmc::platform::cpu::EnergyMethod::DIRECT)
        .value("EWALD", ::pygcmc::platform::cpu::EnergyMethod::EWALD)
        .value("PME", ::pygcmc::platform::cpu::EnergyMethod::PME);
    
    // For now, we skip the GCMCSimulation class bindings since
    // that functionality has been moved to MovementAPI
    // The actual GCMC functionality is exposed through movement module
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc