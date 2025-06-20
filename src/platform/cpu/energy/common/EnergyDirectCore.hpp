#pragma once

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Direct summation for nonbonded interactions
 * 
 * This module provides direct summation implementations for calculating
 * nonbonded interactions between all pairs of atoms.
 */

/**
 * @brief Universal function for calculating all nonbonded interactions
 * 
 * @param state MC state
 * @param use_cutoff Whether to use distance cutoff
 * @param movement_only Whether to calculate only movement residues
 * @param use_pbc Whether to use periodic boundary conditions
 * @param vdw_only Whether to calculate only VDW interactions
 */
void computeNonbondedEnergy(model::MCState& state, bool use_cutoff, bool movement_only, bool use_pbc, bool vdw_only = false);

/**
 * @brief Calculate nonbonded interactions between a single residue and all other active residues
 * 
 * @param state MC state
 * @param residue_idx Index of the residue to calculate
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeResidueNonbondedEnergy(model::MCState& state, int residue_idx, bool use_cutoff, bool use_pbc);

/**
 * @brief Basic direct calculation functions (following existing naming pattern)
 */
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

/**
 * @brief PBC-specific calculation functions
 */
void computeSystemEnergyPBC(model::MCState& state);
void computeSystemEnergyPBCCutoff(model::MCState& state);

/**
 * @brief Unified Direct method interfaces
 */
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 