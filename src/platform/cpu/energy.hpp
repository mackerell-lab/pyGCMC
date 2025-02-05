// src/platform/cpu/energy.hpp

#pragma once

#include "../../model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute a naive (simplified) nonbonded energy based on Lennard-Jones potential,
 *        considering only active movement molecules and all other active molecules.
 *
 * This implementation uses a simplified approach: 
 * - It assumes the distance r = 1.0 for all pairs (i.e. no distance dependence).
 * - It neglects periodic boundary conditions.
 * - For a pair of residues, the Lennard-Jones interaction is computed using the forcefield parameters:
 *      V = eps * [ (Rmin / r)^12  - 2 * (Rmin / r)^6 ]
 *   where:
 *      - eps is obtained from ljEps using index = moveType * maxTypes + residueType
 *      - For simplicity, we assume Rmin = sigma
 * 
 * @param state The current MC state containing residues and force field parameters.
 * @return The computed nonbonded energy.
 */
float computeNaiveNonbondedEnergy(const model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc