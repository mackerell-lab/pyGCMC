// src/platform/cpu/energy.hpp

#pragma once

#include "../../model/montecarlo.hpp"
#include "../platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate nonbonded energies for movement residues only
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between movement residues and all other active residues without distance cutoff.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeMovementResiduesEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues without distance cutoff.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeFullSystemEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with distance cutoff
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues within the specified cutoff distance.
 * 
 * @param state System state containing residues and force field parameters
 */
void computeFullSystemCutoffEnergy(model::MCState& state);

/**
 * @brief Calculate nonbonded energies for the full system with distance cutoff and periodic boundary conditions
 * 
 * This function calculates nonbonded interactions (VDW and electrostatic)
 * between all active residues within the specified cutoff distance,
 * applying periodic boundary conditions using the minimum image convention.
 * 
 * @param state System state containing residues and force field parameters
 * @throws std::runtime_error if box dimensions are invalid for PBC calculation
 */
void computeFullSystemCutoffPBCEnergy(model::MCState& state);

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable);

} // namespace cpu
} // namespace platform
} // namespace pygcmc