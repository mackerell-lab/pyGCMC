// src/platform/cpu/energyDirect.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
extern const float COULOMB;
extern const float MIN_SAFE_DISTANCE;
extern const float MAX_SAFE_ENERGY;

/**
 * @brief Calculate nonbonded energies for movement residues only
 */
void computeMovementEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for movement residues only with distance cutoff
 */
void computeMovementEnergyCutoff(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system
 */
void computeSystemEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with distance cutoff
 */
void computeSystemEnergyCutoff(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with periodic boundary conditions
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues without distance cutoff,
 * applying periodic boundary conditions using the minimum image convention.
 * 
 * @param state System state containing residues and force field parameters
 * @throws std::runtime_error if box dimensions are invalid for PBC calculation
 */
void computeSystemEnergyPBC(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues within the specified cutoff distance,
 * applying periodic boundary conditions using the minimum image convention.
 * 
 * @param state System state containing residues and force field parameters
 * @throws std::runtime_error if box dimensions are invalid for PBC calculation
 */
void computeSystemEnergyPBCCutoff(model::MCState& state);

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 