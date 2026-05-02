#pragma once

#include "model/ModelModule.hpp"
#include "NeighborList.hpp"

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
 * @param vdw_only Whether to calculate only VDW interactions
 */
enum class ResiduePartnerFilter {
    All,
    FixedOnly,
    NonFixedOnly,
};

void computeResidueNonbondedEnergy(model::MCState& state,
                                  int residue_idx,
                                  bool use_cutoff,
                                  bool use_pbc,
                                  bool vdw_only = false,
                                  bool include_pairtypes14_intra = false,
                                  ResiduePartnerFilter partner_filter = ResiduePartnerFilter::All);

/**
 * @brief Compute nonbonded energy (vdw+elec) for a single residue vs the rest, with cutoff and PBC
 *
 * @param state MC state
 * @param residue_idx Index of the residue to calculate
 */
void computeResidueEnergyCutoffPBC(model::MCState& state, int residue_idx);

/**
 * @brief Basic direct calculation functions (following existing naming pattern)
 */
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

/**
 * @brief Movement-specific VDW calculation function
 */
void computeMovementVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

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

/**
 * @brief Compute nonbonded energy using neighbor list
 *
 * @param state MC state
 * @param neighborList Pre-built neighbor list
 * @param use_pbc Whether to use periodic boundary conditions
 * @param vdw_only Whether to calculate only VDW interactions
 */
void computeNonbondedEnergyWithNeighborList(
    model::MCState& state,
    const NeighborList& neighborList,
    bool use_pbc,
    bool vdw_only = false
);

/**
 * @brief Compute residue energy using neighbor list
 *
 * @param state MC state
 * @param residue_idx Index of the residue to calculate
 * @param neighborList Pre-built neighbor list
 * @param use_pbc Whether to use periodic boundary conditions
 * @param vdw_only Whether to calculate only VDW interactions
 */
void computeResidueEnergyWithNeighborList(
    model::MCState& state,
    int residue_idx,
    const NeighborList& neighborList,
    bool use_pbc,
    bool vdw_only = false
);

} // namespace cpu
} // namespace platform
} // namespace pygcmc
