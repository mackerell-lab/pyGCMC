// src/platform/cpu/energy/pgp/PGPPrecompute.hpp

#pragma once

#include "PGPCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Precomputation functions for PGP-PME
 *
 * This module handles parameter setup and grid potential precomputation
 * for the Precomputed Grid-Potential Particle Mesh Ewald algorithm.
 */

// Parameter setup functions
void setPGPParametersImpl(double alpha, const int meshSize[3], double potential_cutoff,
                          const int potentialGridSize[3], int splineOrder, double tolerance);

// Grid potential precomputation functions
void precomputeGridPotentialImpl(model::MCState& state, bool fixed_only = true);

// <agent-hook:pgp_precompute>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
