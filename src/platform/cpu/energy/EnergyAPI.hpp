#pragma once

/**
 * @file EnergyAPI.hpp
 * @brief Public API for energy calculations
 * 
 * This file provides the public interface for energy calculations,
 * extracted from the simulation module to provide direct access to
 * CPU energy computation functions.
 */

#include "../../../model/ModelModule.hpp"
#include "../../../system/common/SystemLogger.hpp"
#include "EnergyModule.hpp"
#include "pme/PMEComposite.hpp"
#include "pme/PMESetup.hpp"
#include "pme/PMEGlobal.hpp"
#include "pgp/PGPCore.hpp"
#include "pgp/PGPComplete.hpp"
#include "pgp/PGPComposite.hpp"
#include "ewald/EwaldComposite.hpp"
#include <cmath>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace energy {

using namespace pygcmc::system::common;
using namespace pygcmc::model::montecarlo;

/**
 * @brief Compute nonbonded energy for movement residues only
 */
inline void computeMovementEnergy(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing nonbonded energy for movement residues");
    }
    platform::cpu::computeMovementEnergy(state, EnergyMethod::DIRECT, false, false);
}

/**
 * @brief Compute nonbonded energy for movement residues with cutoff
 */
inline void computeMovementEnergyCutoff(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing nonbonded energy for movement residues with cutoff");
    }
    platform::cpu::computeMovementEnergyCutoff(state);
}

/**
 * @brief Compute nonbonded energy for the full system
 */
inline void computeSystemEnergy(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing nonbonded energy for all active residues");
    }
    
    platform::cpu::computeSystemEnergy(state, EnergyMethod::DIRECT, false, false);
    
    // Log total energy in debug mode
    if (SystemLogger::isDebugEnabled()) {
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active) {
                total_vdw += state.residues[i].energy_vdw;
                total_elec += state.residues[i].energy_elec;
            }
        }
        
        total_vdw /= 2.0f;
        total_elec /= 2.0f;
        
        SystemLogger::debug("Total system energy: vdw=", total_vdw, 
            ", elec=", total_elec, 
            ", total=", (total_vdw + total_elec));
    }
}

/**
 * @brief Get total energy components from state
 */
inline std::pair<double, double> getTotalEnergyComponents(const MCState& state) {
    double total_elec = 0.0;
    double total_vdw = 0.0;
    
    // Sum energy components from all active residues
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (state.residues[i].active) {
            total_elec += state.residues[i].energy_elec;
            total_vdw += state.residues[i].energy_vdw;
        }
    }
    
    // Divide by 2 to account for double-counting in pairwise interactions
    total_elec /= 2.0;
    total_vdw /= 2.0;
    
    return std::make_pair(total_elec, total_vdw);
}

/**
 * @brief Compute and return total system energy
 */
inline double computeSystemEnergyTotal(MCState& state) {
    computeSystemEnergy(state);
    auto [elec, vdw] = getTotalEnergyComponents(state);
    return elec + vdw;
}

/**
 * @brief Compute nonbonded energy for the full system with cutoff
 */
inline void computeSystemEnergyCutoff(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing cutoff nonbonded energy for all active residues");
    }
    
    platform::cpu::computeSystemEnergyCutoff(state);
    
    // Log total energy in debug mode
    if (SystemLogger::isDebugEnabled()) {
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active) {
                total_vdw += state.residues[i].energy_vdw;
                total_elec += state.residues[i].energy_elec;
            }
        }
        
        total_vdw /= 2.0f;
        total_elec /= 2.0f;
        
        SystemLogger::debug("Total system energy: vdw=", total_vdw, 
            ", elec=", total_elec, 
            ", total=", (total_vdw + total_elec));
    }
}

/**
 * @brief Compute and return total system energy with cutoff
 */
inline double computeSystemEnergyCutoffTotal(MCState& state) {
    platform::cpu::computeSystemEnergyCutoff(state);
    auto [elec, vdw] = getTotalEnergyComponents(state);
    return elec + vdw;
}

/**
 * @brief Compute system energy with periodic boundary conditions
 */
inline void computeSystemEnergyPBC(MCState& state) {
    // Validate box dimensions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing PBC nonbonded energy for all active residues");
    }
    
    platform::cpu::computeSystemEnergyPBC(state);
}

/**
 * @brief Compute system energy with PBC and cutoff
 */
inline void computeSystemEnergyPBCCutoff(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing PBC cutoff nonbonded energy");
    }
    platform::cpu::computeSystemEnergyPBCCutoff(state);
}

/**
 * @brief Compute VdW energy only with cutoff
 */
inline void computeSystemVdwEnergyCutoff(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing VdW energy with cutoff");
    }
    platform::cpu::computeSystemVdwEnergyCutoff(state);
}

// ============================================================================
// Ewald Method Functions
// ============================================================================

/**
 * @brief Set Ewald parameters
 */
inline void setEwaldParameters(float alpha, const int kmax[3], float tolerance = 1e-5f) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Setting Ewald parameters: alpha=", alpha,
            ", kmax=[", kmax[0], ",", kmax[1], ",", kmax[2], "], tolerance=", tolerance);
    }
    
    // Call the actual Ewald implementation
    platform::cpu::setEwaldParameters(alpha, kmax, tolerance);
}

/**
 * @brief Initialize Ewald parameters with automatic optimization
 */
inline void initializeEwaldParameters(float cutoff, const float box[3], 
                                      float alpha = 0.0f, float tolerance = 1e-5f) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Initializing Ewald parameters: cutoff=", cutoff, 
            ", box=[", box[0], ",", box[1], ",", box[2], "]");
    }
    
    // Use EwaldComposite initialization
    double dbox[3] = {box[0], box[1], box[2]};
    EwaldComposite::initialize(cutoff, dbox, alpha, tolerance);
}

/**
 * @brief Compute system energy using Ewald summation
 */
inline void computeSystemEnergyEwald(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using Ewald summation");
    }
    platform::cpu::computeSystemEnergyEwald(state);
}

/**
 * @brief Compute movement energy using Ewald summation
 */
inline void computeMovementEnergyEwald(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing movement energy using Ewald summation");
    }
    platform::cpu::computeMovementEnergyEwald(state);
}

// ============================================================================
// PME Method Functions
// ============================================================================

/**
 * @brief Set PME parameters
 */
inline void setPMEParameters(float alpha, const int meshSize[3], int splineOrder = 4, float tolerance = 1e-5f) {
    // Use PMESetup function - note the namespace is different
    platform::cpu::setPMEParameters(alpha, meshSize, splineOrder, tolerance);
}

/**
 * @brief Initialize PME parameters with automatic optimization
 */
inline void initializePMEParameters(float cutoff, const float box[3], 
                                    float alpha = 0.0f, const int* meshSize = nullptr,
                                    int splineOrder = 4, float tolerance = 1e-5f) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Initializing PME parameters");
    }
    
    // Use PMEComposite initialization
    double dbox[3] = {box[0], box[1], box[2]};
    PMEComposite::initialize(cutoff, dbox, tolerance, alpha, meshSize, splineOrder);
}

/**
 * @brief Compute system energy using PME
 */
inline void computeSystemEnergyPME(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using PME");
    }
    platform::cpu::computeSystemEnergyPME(state);
}

/**
 * @brief Compute movement energy using PME
 */
inline void computeMovementEnergyPME(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing movement energy using PME");
    }
    platform::cpu::computeMovementEnergyPME(state);
}

/**
 * @brief Clear PME engine state
 */
inline void clearPMEEngine() {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Clearing PME engine state");
    }
    platform::cpu::clearPMEState();
}

/**
 * @brief Compute system energy using PME Complete
 */
inline void computeSystemEnergyPMEComplete(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using PME Complete");
    }
    platform::cpu::computeSystemEnergyPMEComplete(state);
}

/**
 * @brief Compute system energy using cutoff Complete
 */
inline void computeSystemEnergyCutoffComplete(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using cutoff Complete");
    }
    platform::cpu::computeSystemEnergyCutoffComplete(state);
}

// ============================================================================
// PGP Method Functions  
// ============================================================================

/**
 * @brief Reset PGP global state
 */
inline void resetPGPState() {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Resetting PGP state");
    }
    // Call platform reset function
    ::pygcmc::platform::cpu::resetPGPState();
}

/**
 * @brief Set PGP parameters
 */
inline void setPGPParameters(float alpha, const int meshSize[3], float potential_cutoff, 
                             const int pairGridSize[3], int splineOrder = 4, float tolerance = 1e-5f) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Setting PGP parameters");
    }
    
    // Call the actual PGP implementation
    ::pygcmc::platform::cpu::setPGPParameters(alpha, meshSize, potential_cutoff, 
                                               pairGridSize, splineOrder, tolerance);
}

/**
 * @brief Initialize PGP parameters with automatic optimization
 */
inline void initializePGPParameters(float cutoff, float pair_cutoff, const float box[3], 
                                    float alpha = 0.0f, const int* meshSize = nullptr,
                                    const int* pairGridSize = nullptr,
                                    int splineOrder = 4, float tolerance = 1e-5f) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Initializing PGP parameters");
    }
    
    // Use default values if not provided
    int default_mesh[3] = {32, 32, 32};
    int default_pair[3] = {64, 64, 64};
    
    const int* mesh = meshSize ? meshSize : default_mesh;
    const int* pair = pairGridSize ? pairGridSize : default_pair;
    
    // Convert float box to double for PGPComposite
    double dbox[3] = {box[0], box[1], box[2]};
    
    // Initialize using PGPComposite
    PGPComposite::initialize(cutoff, dbox, alpha, mesh, pair_cutoff, pair, splineOrder, tolerance);
}

/**
 * @brief Precompute grid potential for PGP
 */
inline void precomputeGridPotential(MCState& state, bool fixed_only) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Precomputing grid potential");
    }
    // Directly call PGPComposite::precomputeGrids to avoid namespace issues
    ::pygcmc::platform::cpu::PGPComposite::precomputeGrids(state, fixed_only);
}

/**
 * @brief Interpolate molecule energy from grid
 */
inline void interpolateMoleculeEnergy(MCState& state, double& energy) {
    ::pygcmc::platform::cpu::interpolateMoleculeEnergy(state, energy);
}

/**
 * @brief Calculate molecule energy from grid
 */
inline double calculateMoleculeEnergy(MCState& state) {
    double energy;
    ::pygcmc::platform::cpu::interpolateMoleculeEnergy(state, energy);
    return energy;
}

/**
 * @brief Compute system energy using PGP
 */
inline void computeSystemEnergyPGP(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using PGP");
    }
    ::pygcmc::platform::cpu::computeSystemEnergyPGP(state);
}

/**
 * @brief Compute movement energy using PGP
 */
inline void computeMovementEnergyPGP(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing movement energy using PGP");
    }
    ::pygcmc::platform::cpu::computeMovementEnergyPGP(state);
}

/**
 * @brief Compute system energy using PGP Fixed
 */
inline void computeSystemEnergyPGPFixed(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using PGP Fixed");
    }
    ::pygcmc::platform::cpu::computeSystemEnergyPGPFixed(state);
}

/**
 * @brief Compute movement energy using PGP Fixed
 */
inline void computeMovementEnergyPGPFixed(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing movement energy using PGP Fixed");
    }
    ::pygcmc::platform::cpu::computeMovementEnergyPGPFixed(state);
}

/**
 * @brief Compute system energy using PGP Complete
 */
inline void computeSystemEnergyPGPComplete(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing system energy using PGP Complete");
    }
    ::pygcmc::platform::cpu::computeSystemEnergyPGPComplete(state);
}

/**
 * @brief Compute movement energy using PGP Complete
 */
inline void computeMovementEnergyPGPComplete(MCState& state) {
    if (SystemLogger::isDebugEnabled()) {
        SystemLogger::debug("Computing movement energy using PGP Complete");
    }
    ::pygcmc::platform::cpu::computeMovementEnergyPGPComplete(state);
}

/**
 * @brief Set energy debug output
 */
inline void setEnergyDebugOutput(bool enable) {
    platform::set_debug_mode(enable);
    if (enable) {
        SystemLogger::setVerbose(true);
        SystemLogger::setLogLevel(system::common::LogLevel::DEBUG);
    }
}

} // namespace energy
} // namespace cpu
} // namespace platform
} // namespace pygcmc