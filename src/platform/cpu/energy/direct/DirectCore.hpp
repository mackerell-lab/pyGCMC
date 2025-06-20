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
 * @brief Universal function for calculating all nonbonded interactions
 * 
 * Core function that handles the direct calculation of non-bonded interactions.
 * This is the main computational engine for direct energy calculations.
 * 
 * @param state System state containing atoms, residues, and forcefield
 * @param use_cutoff Whether to use distance cutoff
 * @param movement_only Whether to calculate only for movement residues
 * @param use_pbc Whether to use periodic boundary conditions
 * @param vdw_only Whether to calculate only VDW interactions
 */
void computeNonbondedEnergy(model::MCState& state, 
                           bool use_cutoff, 
                           bool movement_only = false, 
                           bool use_pbc = false, 
                           bool vdw_only = false);

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 