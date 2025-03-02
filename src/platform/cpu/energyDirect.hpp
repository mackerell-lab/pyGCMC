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
 * @brief Calculate system non-bonded energy using direct calculation method
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

/**
 * @brief Calculate non-bonded energy for movement residues using direct calculation method
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

/**
 * @brief Calculate only van der Waals energy for the system (with cutoff)
 * 
 * This function only calculates van der Waals interactions, not electrostatic interactions.
 * Mainly used in conjunction with the Ewald summation method, where electrostatic interactions are handled separately.
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff (generally true)
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

// Functions retained for compatibility with old interfaces
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemEnergyPBC(model::MCState& state);
void computeSystemEnergyPBCCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 