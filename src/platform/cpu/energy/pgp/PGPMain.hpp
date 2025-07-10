#pragma once

/**
 * @brief PGP Module Unified Entry Point - Precomputed Grid-Potential Module
 * 
 * This file aggregates all functionality of the PGP module, external code only needs to include this file.
 * 
 * Functional components:
 * - PGPCore: Core parameters and initialization
 * - PGPGrid: Grid operations and potential grids
 * - PGPInterpolation: Grid interpolation and energy calculation
 * - PGPPrecompute: Precomputation algorithm optimization
 * - PGPReal: Real-space energy calculation
 * - PGPSelf: Self-energy correction
 * - PGPSystem: Complete system energy evaluation
 * - PGPComposite: Advanced unified interface
 * 
 * Typical usage:
 *   #include "pgp/PGPMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyPGP(state);
 *   computeMovementEnergyPGP(state);
 * 
 * @note Module functionality guide:
 * - Energy calculation: PGPSystem.hpp -> computeSystemEnergyPGP, computeMovementEnergyPGP
 * - Grid interpolation: PGPSystem.hpp -> interpolateMoleculeEnergy, calculateMoleculeEnergy
 * - Grid precomputation: PGPPrecompute.hpp -> precomputeGridPotential, setPGPParameters
 * - Real space: PGPReal.hpp -> computeRealSpacePGP
 * - Self energy: PGPSelf.hpp -> computeSelfEnergyPGP
 * - Grid operations: PGPGrid.hpp -> initializePotentialGrid
 * - Parameter setting: PGPCore.hpp -> setPGPParameters, PGPParams structure
 */

// Aggregate all sub-functions of the PGP module
#include "PGPCore.hpp"
#include "PGPGrid.hpp"
#include "PGPInterpolation.hpp"
#include "PGPPrecompute.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "PGPSystem.hpp"
#include "PGPComposite.hpp"
#include "PGPComplete.hpp"

// <agent-hook:pgp_main>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PGP Module Information
 */
namespace PGPInfo {
    constexpr const char* VERSION = "1.0.0";
    constexpr const char* DESCRIPTION = "Modular Precomputed Grid-Potential Implementation";
    constexpr int NUM_MODULES = 7;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "PGPCore",
        "PGPGrid",
        "PGPInterpolation", 
        "PGPPrecompute",
        "PGPRealSpace",
        "PGPSelfEnergy",
        "PGPSystemEnergy"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 