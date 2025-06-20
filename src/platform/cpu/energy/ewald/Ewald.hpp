#pragma once

/**
 * @brief Ewald Summation Module
 * 
 * This header provides a unified interface to the modular Ewald implementation.
 * The Ewald module is split into the following components:
 * 
 * - EwaldCore: Core parameters and initialization
 * - EwaldRealSpace: Real space energy calculations
 * - EwaldReciprocal: Reciprocal space energy calculations
 * - EwaldSelf: Self-energy corrections
 * - EwaldComposite: High-level unified interface
 * 
 * For most users, including EwaldComposite.hpp is sufficient, as it provides
 * backward-compatible interfaces and orchestrates all other modules.
 */

// Core functionality
#include "EwaldCore.hpp"

// Individual modules
#include "EwaldRealSpace.hpp"
#include "EwaldReciprocal.hpp"
#include "EwaldSelf.hpp"

// Unified interface
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