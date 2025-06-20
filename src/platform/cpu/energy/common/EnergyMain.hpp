#pragma once

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