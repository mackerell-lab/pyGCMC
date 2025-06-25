#pragma once

/**
 * @brief Ewald Module Unified Entry Point - Ewald Summation Module
 * 
 * This file aggregates all functionality of the Ewald module, external code only needs to include this file.
 * 
 * Functional components:
 * - EwaldCore: Core parameters and initialization
 * - EwaldReal: Real-space energy calculation
 * - EwaldRecip: Reciprocal-space energy calculation
 * - EwaldSelf: Self-energy correction
 * - EwaldComposite: Advanced unified interface
 * 
 * Typical usage:
 *   #include "ewald/EwaldMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyEwald(state);
 *   computeMovementEnergyEwald(state);
 * 
 * @note Module functionality guide:
 * - Energy calculation: EwaldComposite.hpp -> computeSystemEnergyEwald, computeMovementEnergyEwald
 * - Real space: EwaldReal.hpp -> computeRealSpaceEwald, calcPairEnergyEwaldRealSpace
 * - Reciprocal space: EwaldRecip.hpp -> computeReciprocalEnergy
 * - Self energy: EwaldSelf.hpp -> computeSelfEnergy
 * - Initialization: EwaldInterface.hpp -> initializeEwald, isEwaldInitialized
 * - Parameter setting: EwaldCore.hpp -> setEwaldParameters, autoAdjustParameters
 */

// Aggregate all sub-functions of the Ewald module
#include "EwaldCore.hpp"
#include "EwaldReal.hpp" 
#include "EwaldRecip.hpp"
#include "EwaldSelf.hpp"
#include "EwaldComposite.hpp"

// <agent-hook:ewald_main>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Ewald Module Information
 */
namespace EwaldInfo {
    constexpr const char* VERSION = "1.0.0";
    constexpr const char* DESCRIPTION = "Modular Ewald Summation Implementation";
    constexpr int NUM_MODULES = 4;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "EwaldCore",
        "EwaldRealSpace",
        "EwaldReciprocal", 
        "EwaldSelf"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 