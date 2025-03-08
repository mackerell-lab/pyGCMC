// src/platform/cpu/energy.hpp

#pragma once

// First include common interfaces
#include "energyCommon.hpp"

// Then include implementations
#include "energyDirect.hpp"
#include "energyEwald.hpp"
#include "energyPME.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief CPU Platform Energy Calculation Module
 * 
 * This module provides energy calculation functionality on the CPU platform, including:
 * 1. Direct calculation method (Direct): Uses explicit summation to calculate van der Waals and Coulomb interactions
 * 2. Ewald summation method: Optimized calculation for long-range electrostatic interactions in periodic systems
 * 3. Particle Mesh Ewald (PME) method: Fast approximation of Ewald summation using FFT for larger systems
 * 
 * Basic usage examples:
 * 
 * // Direct calculation (no cutoff, no PBC):
 * computeSystemEnergy(state, EnergyMethod::DIRECT, false, false);
 * 
 * // Direct calculation (with cutoff, with PBC):
 * computeSystemEnergy(state, EnergyMethod::DIRECT, true, true);
 * 
 * // Ewald summation (automatically uses PBC and cutoff):
 * computeSystemEnergy(state, EnergyMethod::EWALD);
 * 
 * // Particle Mesh Ewald (automatically uses PBC and cutoff):
 * computeSystemEnergy(state, EnergyMethod::PME);
 * 
 * Note: When using the Ewald method, Ewald parameters must be initialized first:
 * initializeEwaldParameters(cutoff, box);
 * 
 * Note: When using the PME method, PME parameters must be initialized first:
 * initializePMEParameters(cutoff, box);
 */

/**
 * @brief Get the name of the energy calculation method
 * 
 * @param method Energy calculation method enum
 * @return String name of the method
 */
inline std::string getEnergyMethodName(EnergyMethod method) {
    switch (method) {
        case EnergyMethod::DIRECT: return "Direct";
        case EnergyMethod::EWALD: return "Ewald";
        case EnergyMethod::PME: return "PME";
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
inline double getTotalEnergy(const model::MCState& state, EnergyMethod method) {
    if (method == EnergyMethod::EWALD || method == EnergyMethod::PME) {
        return getEwaldTotalEnergy(state);
    } else {
        double total = 0.0;
        for (const auto& residue : state.residues) {
            if (residue.active) {
                total += residue.energy_vdw + residue.energy_elec;
            }
        }
        return total;
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 