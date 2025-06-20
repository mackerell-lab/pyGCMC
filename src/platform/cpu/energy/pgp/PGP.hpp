#pragma once

/**
 * @brief PGP (Precomputed Grid-Potential) Module
 * 
 * This header provides a unified interface to the modular PGP implementation.
 * The PGP module is split into the following components:
 * 
 * - PGPCore: Core parameters and initialization
 * - PGPGrid: Grid operations and potential grids
 * - PGPPrecompute: Precomputation algorithms for optimization
 * - PGPRealSpace: Real space energy calculations
 * - PGPSelfEnergy: Self energy corrections
 * - PGPSystemEnergy: Complete system energy evaluation and grid interpolation
 * - PGPComposite: High-level unified interface
 * 
 * For most users, including PGPComposite.hpp is sufficient, as it provides
 * backward-compatible interfaces and orchestrates all other modules.
 * 
 * @brief Function Location Guide for AI Agents:
 * - Energy calculation: PGPSystemEnergy.hpp -> computeSystemEnergyPGP, computeMovementEnergyPGP
 * - Grid interpolation: PGPSystemEnergy.hpp -> interpolateMoleculeEnergy, calculateMoleculeEnergy
 * - Grid precomputation: PGPPrecompute.hpp -> precomputeGridPotential, setPGPParameters
 * - Real space: PGPRealSpace.hpp -> computeRealSpacePGP
 * - Self energy: PGPSelfEnergy.hpp -> computeSelfEnergyPGP
 * - Grid operations: PGPGrid.hpp -> initializePotentialGrid
 * - Parameter setup: PGPCore.hpp -> setPGPParameters, PGPParams structure
 */

// Core functionality
#include "PGPCore.hpp"

// Individual modules
#include "PGPGrid.hpp"
#include "PGPPrecompute.hpp"
#include "PGPRealSpace.hpp"
#include "PGPSelfEnergy.hpp"
#include "PGPSystemEnergy.hpp"

// Unified interface
#include "PGPComposite.hpp"

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