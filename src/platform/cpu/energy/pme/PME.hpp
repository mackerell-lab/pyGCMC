#pragma once

/**
 * @brief PME (Particle Mesh Ewald) Module
 * 
 * This header provides a unified interface to the modular PME implementation.
 * The PME module is split into the following components:
 * 
 * - PMECore: Core parameters and initialization
 * - PMEFFT: Custom FFT implementation
 * - PMESpline: B-spline interpolation functions
 * - PMEGrid: Grid operations and charge spreading
 * - PMEReciprocal: Reciprocal space energy calculations
 * - PMERealSpace: Real space energy calculations
 * - PMESelf: Self-energy corrections
 * - PMEComposite: High-level unified interface
 * 
 * For most users, including PMEComposite.hpp is sufficient, as it provides
 * backward-compatible interfaces and orchestrates all other modules.
 */

// Core functionality
#include "PMECore.hpp"

// Individual modules
#include "PMEFFT.hpp"
#include "PMESpline.hpp"
#include "PMEGrid.hpp"
#include "PMEReciprocal.hpp"
#include "PMERealSpace.hpp"
#include "PMESelf.hpp"

// Unified interface
#include "PMEComposite.hpp"

// <agent-hook:pme_main>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PME Module Information
 */
namespace PMEInfo {
    constexpr const char* VERSION = "1.0.0";
    constexpr const char* DESCRIPTION = "Modular Particle Mesh Ewald Implementation";
    constexpr int NUM_MODULES = 7;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "PMECore",
        "PMEFFT", 
        "PMESpline",
        "PMEGrid",
        "PMEReciprocal",
        "PMERealSpace",
        "PMESelf"
    };
}

} // namespace cpu
} // namespace platform  
} // namespace pygcmc 