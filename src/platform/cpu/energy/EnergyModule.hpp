#pragma once

#include <string>

// Include all energy calculation modules
#include "common/EnergySystemInterface.hpp"
#include "common/DirectSummation.hpp"
#include "lj/LJSwitching.hpp"
#include "coulomb/PairEnergyCalculation.hpp"
#include "ewald/EwaldComposite.hpp"
#include "pme/PMEComposite.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace energy {

/**
 * @brief CPU Platform Energy Calculation Module
 * 
 * This module provides energy calculation functionality on the CPU platform, including:
 * 1. Direct calculation method (Direct): Uses explicit summation to calculate van der Waals and Coulomb interactions
 * 2. Ewald summation method: Optimized calculation for long-range electrostatic interactions in periodic systems
 * 3. Particle Mesh Ewald (PME) method: Fast approximation of Ewald summation using FFT for larger systems
 * 
 * Modern usage examples:
 * 
 * // Direct calculation using the new modular interface:
 * computeSystemEnergyDirect(state, false, false);
 * 
 * // Direct calculation with cutoff and PBC:
 * computeSystemEnergyDirect(state, true, true);
 * 
 * // Using the unified interface:
 * computeSystemEnergy(state, EnergyMethod::DIRECT, true, true);
 * computeSystemEnergy(state, EnergyMethod::EWALD);
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
        return pygcmc::platform::cpu::getTotalEnergy(state);
    }
}

} // namespace energy
} // namespace cpu
} // namespace platform
} // namespace pygcmc 