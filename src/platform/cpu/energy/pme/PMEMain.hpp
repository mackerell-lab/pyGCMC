#pragma once

/**
 * @brief PME Module Unified Entry Point - Particle Mesh Ewald Module
 * 
 * This file aggregates all functionality of the PME module, external code only needs to include this file.
 * 
 * Functional components:
 * - PMECore: Core parameter structures and basic functions
 * - PMEConfig: Parameter setting and configuration management
 * - PMESetup: Initialization functions
 * - PMEInterface: Advanced energy calculation interface
 * - PMEFFT: Custom FFT implementation
 * - PMESpline: B-spline interpolation functions
 * - PMEGrid: Grid operations and charge distribution
 * - PMEReal: Real-space energy calculation
 * - PMERecip: Reciprocal-space energy calculation
 * - PMESelf: Self-energy correction
 * - PMEComposite: Unified interface and convenience functions
 * 
 * Typical usage:
 *   #include "pme/PMEMain.hpp"
 *   
 *   using namespace pygcmc::platform::cpu;
 *   computeSystemEnergyPME(state);
 *   computeMovementEnergyPME(state);
 * 
 * @note Module functionality guide:
 * - Energy calculation: PMEInterface.hpp -> computeSystemEnergyPME, computeMovementEnergyPME
 * - Component calculation: PMEInterface.hpp -> computeReciprocalPME, computeSelfEnergyPME, computeRealSpacePME
 * - Parameter setting: PMESetup.hpp -> setPMEParameters, autoAdjustPMEParameters
 * - Initialization: PMESetup.hpp -> initializePMEParameters, initializePMETables, initializePMEBsplines
 * - Real space: PMEReal.hpp -> computeRealSpaceEnergy, calcPairEnergyPME
 * - Self energy: PMESelf.hpp -> computeSelfEnergyPME, calculateParticleSelfEnergy
 * - Reciprocal space: PMERecip.hpp -> reciprocal space calculation
 * - Grid operations: PMEGrid.hpp, PMEGridMap.hpp -> grid management and charge distribution
 * - Parameter structures: PMECore.hpp -> PMEParams structure and basic functions
 */

// Aggregate all sub-functions of the PME module
#include "PMECore.hpp"
#include "PMEConfig.hpp"
#include "PMESetup.hpp"
#include "PMEInterface.hpp"
#include "PMESpline.hpp"
#include "PMEGrid.hpp"
#include "PMEReal.hpp"
#include "PMERecip.hpp"
#include "PMESelf.hpp"
#include "PMEComposite.hpp"
#include "PMEComplete.hpp"

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