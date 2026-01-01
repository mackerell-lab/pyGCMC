// src/bindings/platform/PlatformEnergyBindings.cpp
// Python bindings for platform energy calculation functions

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "platform/cpu/energy/common/EnergyDirectCore.hpp"
#include "model/montecarlo/MCMain.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace platform {

void init_energy_bindings(py::module& m) {
    // Create platform and cpu submodules if they don't exist
    py::module m_platform;
    py::module m_cpu;
    
    // Check if platform submodule already exists
    if (!py::hasattr(m, "platform")) {
        m_platform = m.def_submodule("platform", "Platform-specific implementations");
    } else {
        m_platform = m.attr("platform");
    }
    
    // Check if cpu submodule already exists
    if (!py::hasattr(m_platform, "cpu")) {
        m_cpu = m_platform.def_submodule("cpu", "CPU implementations");
    } else {
        m_cpu = m_platform.attr("cpu");
    }
    
    // Expose energy calculation functions
    m_cpu.def("computeSystemEnergyPBCCutoff",
              &pygcmc::platform::cpu::computeSystemEnergyPBCCutoff,
              py::arg("state"),
              "Compute full system energy with PBC and cutoff");
    
    m_cpu.def("computeResidueEnergyCutoffPBC",
              &pygcmc::platform::cpu::computeResidueEnergyCutoffPBC,
              py::arg("state"),
              py::arg("residue_idx"),
              "Compute single residue interaction energy with PBC and cutoff");
    
    // Also expose other useful energy functions
    m_cpu.def("computeSystemEnergy",
              &pygcmc::platform::cpu::computeSystemEnergy,
              py::arg("state"),
              "Compute full system energy without cutoff");
    
    m_cpu.def("computeSystemEnergyCutoff",
              &pygcmc::platform::cpu::computeSystemEnergyCutoff,
              py::arg("state"),
              "Compute full system energy with cutoff but no PBC");
    
    m_cpu.def("computeSystemEnergyPBC",
              &pygcmc::platform::cpu::computeSystemEnergyPBC,
              py::arg("state"),
              "Compute full system energy with PBC but no cutoff");
    
    m_cpu.def("computeMovementEnergyCutoff",
              &pygcmc::platform::cpu::computeMovementEnergyCutoff,
              py::arg("state"),
              "Compute movement residues energy with cutoff");

    py::enum_<pygcmc::platform::cpu::ResiduePartnerFilter>(m_cpu, "ResiduePartnerFilter")
        .value("All", pygcmc::platform::cpu::ResiduePartnerFilter::All)
        .value("FixedOnly", pygcmc::platform::cpu::ResiduePartnerFilter::FixedOnly)
        .value("NonFixedOnly", pygcmc::platform::cpu::ResiduePartnerFilter::NonFixedOnly)
        .export_values();
    
    m_cpu.def("computeResidueNonbondedEnergy",
              &pygcmc::platform::cpu::computeResidueNonbondedEnergy,
              py::arg("state"),
              py::arg("residue_idx"),
              py::arg("use_cutoff"),
              py::arg("use_pbc"),
              py::arg("vdw_only") = false,
              py::arg("include_pairtypes14_intra") = false,
              py::arg("partner_filter") = pygcmc::platform::cpu::ResiduePartnerFilter::All,
              "Compute residue nonbonded energy with flexible options");
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc
