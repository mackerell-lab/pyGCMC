#pragma once

/**
 * @brief PME (Particle Mesh Ewald) Module
 * 
 * This header provides a unified interface to the modular PME implementation.
 * The PME module is split into the following components:
 * 
 * - PMECore: Core parameter structures and basic functions
 * - PMESetup: Parameter setting and initialization functions
 * - PMEInterface: High-level energy calculation interfaces
 * - PMEFFT: Custom FFT implementation
 * - PMESpline: B-spline interpolation functions
 * - PMEGrid: Grid operations and charge spreading
 * - PMEReciprocal: Reciprocal space energy calculations
 * - PMERealSpace: Real space energy calculations
 * - PMESelf: Self-energy corrections
 * - PMEComposite: Unified interface and convenience functions
 * 
 * For most users, including PMEComposite.hpp is sufficient, as it provides
 * backward-compatible interfaces and orchestrates all other modules.
 * 
 * @brief Function Location Guide for AI Agents:
 * - Energy calculation: PMEInterface.hpp -> computeSystemEnergyPME, computeMovementEnergyPME
 * - Component calculations: PMEInterface.hpp -> computeReciprocalPME, computeSelfEnergyPME, computeRealSpacePME
 * - Parameter setting: PMESetup.hpp -> setPMEParameters, autoAdjustPMEParameters
 * - Initialization: PMESetup.hpp -> initializePMEParameters, initializePMETables, initializePMEBsplines
 * - Real space: PMERealSpace.hpp -> computeRealSpaceEnergy, calcPairEnergyPME
 * - Self energy: PMESelf.hpp -> computeSelfEnergyPME, calculateParticleSelfEnergy
 * - Reciprocal space: PMEReciprocal.hpp -> reciprocal space calculations
 * - Grid operations: PMEGrid.hpp, PMEGridMapping.hpp -> grid management and charge spreading
 * - Parameter structure: PMECore.hpp -> PMEParams structure and basic functions
 */

// Core functionality
#include "PMECore.hpp"
#include "PMESetup.hpp"
#include "PMEInterface.hpp"

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
    constexpr int NUM_MODULES = 9;
    
    // Module names for debugging and introspection
    constexpr const char* MODULE_NAMES[] = {
        "PMECore",
        "PMESetup",
        "PMEInterface",
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