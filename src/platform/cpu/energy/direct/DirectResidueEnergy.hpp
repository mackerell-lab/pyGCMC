#pragma once

#include "../common/EnergyConstants.hpp"
#include "../common/EnergyUtils.hpp"
#include "../common/EnergyPairCalculation.hpp"
#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

/**
 * @brief Calculate nonbonded interactions between a single residue and all other active residues
 * 
 * This function calculates the energy contribution of one residue interacting with all other
 * active residues in the system. It's the building block for both full system and movement
 * residue energy calculations.
 * 
 * @param state System state containing atoms, residues, and forcefield
 * @param residue_idx Index of the residue to calculate energy for
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeResidueNonbondedEnergy(
    model::MCState& state,
    int residue_idx,
    bool use_cutoff = false,
    bool use_pbc = false
);

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 