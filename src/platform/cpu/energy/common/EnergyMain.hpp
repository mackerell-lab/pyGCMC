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

} // namespace energy
} // namespace cpu
} // namespace platform
} // namespace pygcmc 