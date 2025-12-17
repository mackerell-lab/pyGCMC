#pragma once

#include <string>

/**
 * @brief Energy Module Unified Entry Point - Common Energy Calculation Module
 * 
 * This file aggregates all functionality of the common energy calculation module, external code only needs to include this file.
 * 
 * Features include:
 * 1. Unified energy calculation interface (EnergyInterface)
 * 2. Direct summation for nonbonded interactions (EnergyDirectCore)
 * 3. Support for cutoff, PBC, and various calculation modes
 * 4. Energy method enumeration and unified interfaces
 * 
 * Typical usage:
 *   #include "common/EnergyMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyDirect(state, true, true);
 *   computeSystemEnergy(state, EnergyMethod::DIRECT, true, false);
 * 
 * @note This is the unified interface of the common energy calculation module, external modules should include this instead of individual header files
 */

// Aggregate all sub-functions of the common energy calculation module
#include "EnergyInterface.hpp"
#include "EnergyDirectCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace energy {

/**
 * @brief Get the name of the energy calculation method
 * 
 * @param method Energy calculation method enum
 * @return String name of the method
 */
inline std::string getEnergyMethodName(pygcmc::platform::cpu::EnergyMethod method) {
    switch (method) {
        case pygcmc::platform::cpu::EnergyMethod::DIRECT: return "Direct";
        case pygcmc::platform::cpu::EnergyMethod::EWALD: return "Ewald";
        case pygcmc::platform::cpu::EnergyMethod::PME: return "PME";
        default: return "Unknown";
    }
}

/**
 * @brief Get the total energy of the current system
 * 
 * @param state System state
 * @param method Energy calculation method to use
 * @return Total energy (kJ/mol)
 */
inline double getTotalEnergy(const model::MCState& state, pygcmc::platform::cpu::EnergyMethod method) {
    if (method == pygcmc::platform::cpu::EnergyMethod::EWALD || method == pygcmc::platform::cpu::EnergyMethod::PME) {
        return pygcmc::platform::cpu::getEwaldTotalEnergy(state);
    } else {
        return pygcmc::platform::cpu::getTotalEnergy(state);
    }
}

/**
 * @brief Get total energy using a unique-pairs convention for pair interactions
 *
 * Semantics:
 * - DIRECT: residue energies are stored one-to-all (double-counted for pair terms),
 *   so total unique energy is 0.5 * sum(residue energies).
 * - EWALD/PME: electrostatics total is stored in state.ewald_energy.total (real+reciprocal+self),
 *   while VdW residue energies are one-to-all; total unique energy is
 *   state.ewald_energy.total + 0.5 * sum(VdW residue energies).
 *
 * This keeps the existing double-counted APIs intact while providing a safe total energy
 * for ΔU = E_after - E_before paths (e.g., GCMC insertion/deletion in EWALD/PME modes).
 */
inline double getTotalEnergyUniquePairs(const model::MCState& state, pygcmc::platform::cpu::EnergyMethod method) {
    double vdwSum = 0.0;
    double elecSum = 0.0;
    for (const auto& residue : state.residues) {
        if (!residue.active) continue;
        vdwSum += residue.energy_vdw;
        elecSum += residue.energy_elec;
    }

    if (method == pygcmc::platform::cpu::EnergyMethod::EWALD || method == pygcmc::platform::cpu::EnergyMethod::PME) {
        // Electrostatics total is tracked separately for Ewald/PME; residue electrostatics
        // are only a partitioning of the real-space part and should not be halved here.
        return state.ewald_energy.total + 0.5 * vdwSum;
    }

    // DIRECT: both VdW and electrostatics are stored as one-to-all (double-counted).
    return 0.5 * (vdwSum + elecSum);
}

} // namespace energy
} // namespace cpu
} // namespace platform
} // namespace pygcmc
