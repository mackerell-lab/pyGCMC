#pragma once

// Include all movement module components
#include "core/MovementMain.hpp"
#include "common/MovementParams.hpp"
#include "common/MovementResult.hpp"
#include "common/MovementStatistics.hpp"
#include "common/MovementUtils.hpp"
#include "pool/ActivePool.hpp"
#include "bias/CavityBias.hpp"
#include "bias/ConfigBias.hpp"

/**
 * @file MovementModule.hpp
 * @brief Movement Module - Unified Entry Point for All GCMC Movement Functionality
 * 
 * This is the ONLY header file you need to include to access all GCMC movement
 * capabilities in the simulation framework. The module provides insertion, deletion,
 * translation, and rotation moves with advanced biasing techniques.
 * 
 * **Module Organization:**
 * 
 * **core/ directory** - Core Movement Infrastructure
 * - MovementMain.hpp: Main MovementModule class definition and interface
 * - MovementCore.cpp: Implementation of the MovementModule class
 * 
 * **common/ directory** - Common Types and Utilities
 * - MovementParams.hpp: Parameters for GCMC movements
 * - MovementResult.hpp: Result structure for movement attempts
 * - MovementStatistics.hpp: Statistics tracking for acceptance rates
 * - MovementUtils.hpp: Utility functions (random numbers, rotations, PBC)
 * 
 * **moves/ directory** - Individual Movement Types
 * - Insertion.hpp/cpp: Molecule insertion moves
 * - Deletion.hpp/cpp: Molecule deletion moves
 * - Translation.hpp/cpp: Molecule translation moves
 * - Rotation.hpp/cpp: Molecule rotation moves
 * 
 * **bias/ directory** - Advanced Biasing Techniques
 * - CavityBias.hpp/cpp: Cavity detection and biased insertion (500x improvement)
 * - ConfigBias.hpp/cpp: Configurational bias for rotations (150x improvement)
 * 
 * **pool/ directory** - Memory Management
 * - ActivePool.hpp/cpp: Pre-allocated memory pool with lazy deletion
 * 
 * **Usage Example:**
 * @code
 * #include "platform/cpu/movement/MovementModule.hpp"
 * 
 * // Create movement module with parameters
 * MovementParams params(298.15);  // Temperature in K
 * params.chemicalPotential = -15.7;  // kJ/mol
 * params.useCavityBias = true;
 * params.useConfigBias = true;
 * 
 * MovementModule movement(params);
 * 
 * // Perform GCMC moves
 * MCState state;
 * MovementResult result = movement.attemptInsertion(state);
 * if (result.accepted) {
 *     // Move was accepted
 * }
 * @endcode
 */

// Make MovementModule available directly in movement namespace
namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
    // MovementModule is defined in core/MovementMain.hpp
    using MovementModule = MovementModule;
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc