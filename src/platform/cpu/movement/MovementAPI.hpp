#pragma once

/**
 * @file MovementAPI.hpp
 * @brief Public API for GCMC movement operations
 * 
 * This file provides the public interface for movement operations,
 * extracted from the simulation module to provide direct access to
 * CPU movement functions.
 */

#include "../../../model/ModelModule.hpp"
#include "../../../system/common/SystemLogger.hpp"
#include "MovementModule.hpp"
#include "gcmc/GCMCModule.hpp"
#include "reservoir/fragment_reservoir.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace pygcmc::system::common;
using namespace pygcmc::model::montecarlo;

// Re-export key classes for convenience
using GCMCModule = pygcmc::platform::cpu::movement::gcmc::GCMCModule;
using GCMCConfig = pygcmc::platform::cpu::movement::gcmc::GCMCModule::Config;
using FragmentReservoir = pygcmc::platform::cpu::movement::FragmentReservoir;

// MovementResult is already defined in MovementModule.hpp

/**
 * @brief Global GCMC module instance (for backward compatibility)
 */
inline GCMCModule* getGlobalGCMCModule() {
    static GCMCModule globalModule;
    return &globalModule;
}

/**
 * @brief Global Fragment Reservoir instance (for backward compatibility)
 */
inline FragmentReservoir* getGlobalFragmentReservoir() {
    static FragmentReservoir globalReservoir;
    return &globalReservoir;
}

// ============================================================================
// Convenience Functions (extracted from simulation)
// ============================================================================

/**
 * @brief Initialize GCMC module with configuration
 */
inline void initializeGCMC(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Initializing GCMC module");
    }
    getGlobalGCMCModule()->initialize(state);
}

/**
 * @brief Perform a single GCMC move
 */
inline bool performGCMCMove() {
    return getGlobalGCMCModule()->performMove();
}

/**
 * @brief Run multiple GCMC steps
 */
inline void runGCMCSteps(int nSteps) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Running ", nSteps, " GCMC steps");
    }
    
    getGlobalGCMCModule()->runSteps(nSteps);
}

/**
 * @brief Get GCMC statistics
 */
inline std::map<std::string, double> getGCMCStatistics() {
    // Convert GCMCStats to map - for now just return empty map
    // auto stats = getGlobalGCMCModule()->getStatistics();
    return std::map<std::string, double>();
}

// ============================================================================
// Fragment Reservoir Functions
// ============================================================================

/**
 * @brief Add a fragment to the reservoir
 */
inline void addFragmentToReservoir(const std::string& name, 
                                   const std::vector<MCAtom>& atoms,
                                   double chemicalPotential = -15.7) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Adding fragment '", name, "' to reservoir with ", 
                           atoms.size(), " atoms and chemical potential ", chemicalPotential);
    }
    
    // Create and add fragment to the reservoir
    // Note: Actual implementation depends on FragmentTemplate structure
    // For now, just acknowledge the parameters to avoid warnings
    (void)name;
    (void)atoms;
    (void)chemicalPotential;
    
    // TODO: Implement when FragmentTemplate is available
    // FragmentTemplate tmpl;
    // tmpl.name = name;
    // tmpl.atoms = atoms;
    // getGlobalFragmentReservoir()->addTemplate(tmpl);
}

/**
 * @brief Clear the fragment reservoir
 */
inline void clearFragmentReservoir() {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Clearing fragment reservoir");
    }
    // For now, just log - actual implementation would clear the reservoir
    // getGlobalFragmentReservoir()->clear();
}

/**
 * @brief Get number of fragments in reservoir
 */
inline size_t getFragmentCount() {
    // For now, return 0 - actual implementation would query the reservoir
    // return getGlobalFragmentReservoir()->size();
    return 0;
}

// ============================================================================
// Direct Movement Functions
// ============================================================================

/**
 * @brief Perform insertion move
 */
inline bool attemptInsertion(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Attempting insertion");
    }
    // Ensure state is initialized in GCMC module
    getGlobalGCMCModule()->initialize(state);
    return getGlobalGCMCModule()->performMove();
}

/**
 * @brief Perform deletion move
 */
inline bool attemptDeletion(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Attempting deletion");
    }
    // Ensure state is initialized in GCMC module
    getGlobalGCMCModule()->initialize(state);
    return getGlobalGCMCModule()->performMove();
}

/**
 * @brief Perform translation move
 */
inline bool attemptTranslation(MCState& state, int residueIndex) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Attempting translation of residue ", residueIndex);
    }
    // Ensure state is initialized in GCMC module
    getGlobalGCMCModule()->initialize(state);
    // TODO: Pass residueIndex to the move selector when available
    (void)residueIndex; // Acknowledge parameter to avoid warning
    return getGlobalGCMCModule()->performMove();
}

/**
 * @brief Perform rotation move
 */
inline bool attemptRotation(MCState& state, int residueIndex) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Attempting rotation of residue ", residueIndex);
    }
    // Ensure state is initialized in GCMC module  
    getGlobalGCMCModule()->initialize(state);
    // TODO: Pass residueIndex to the move selector when available
    (void)residueIndex; // Acknowledge parameter to avoid warning
    return getGlobalGCMCModule()->performMove();
}

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Set random seed for movement operations
 */
inline void setRandomSeed(unsigned int seed) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Setting random seed to ", seed);
    }
    // This would need to be implemented in MovementModule
    // For now, just log the action
}

/**
 * @brief Enable/disable cavity bias globally
 */
inline void setCavityBiasEnabled(bool enabled) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Setting cavity bias ", enabled ? "enabled" : "disabled");
    }
    // This would modify default params
}

/**
 * @brief Enable/disable configurational bias globally
 */
inline void setConfigBiasEnabled(bool enabled) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Setting config bias ", enabled ? "enabled" : "disabled");
    }
    // This would modify default params
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc