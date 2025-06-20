#pragma once

#include "PGPCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Interpolation functions for PGP-PME
 * 
 * This module handles energy interpolation from precomputed potential grids
 * using B-spline interpolation methods.
 */

// Core interpolation functions
void interpolateMoleculeEnergyImpl(model::MCState& state, double& energy);

double calculateMoleculeEnergyImpl(model::MCState& state);

double computeMoleculeEnergyGlobalImpl(model::MCState& state, const std::vector<int>& movementResidues, 
                                       const std::vector<int>& nearbyResidues, int threadIndex);

// <agent-hook:pgp_interpolation>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 