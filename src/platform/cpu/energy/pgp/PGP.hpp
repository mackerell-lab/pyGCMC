#pragma once

/**
 * @brief PGP (Precomputed Grid-Potential) Module
 * 
 * This header provides a unified interface to the modular PGP implementation.
 * The PGP module is split into the following components:
 * 
 * - PGPCore: Core parameters and initialization
 * - PGPGrid: Grid operations and potential grids
 * - PGPInterpolation: Interpolation functions for grid-based calculations
 * - PGPPrecompute: Precomputation algorithms for optimization
 * - PGPEvaluator: Energy evaluation functions
 * - PGPComposite: High-level unified interface
 * 
 * For most users, including PGPComposite.hpp is sufficient, as it provides
 * backward-compatible interfaces and orchestrates all other modules.
 */

// Core functionality
#include "PGPCore.hpp"

// Individual modules
#include "PGPGrid.hpp"
#include "PGPInterpolation.hpp"
#include "PGPPrecompute.hpp"
#include "PGPEvaluator.hpp"

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
    constexpr int NUM_MODULES = 5;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "PGPCore",
        "PGPGrid",
        "PGPInterpolation", 
        "PGPPrecompute",
        "PGPEvaluator"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 