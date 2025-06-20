#pragma once

#include "PGPCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Energy evaluation functions for PGP-PME
 * 
 * This module handles various energy calculations including real space,
 * self energy correction, and complete system energy evaluation.
 */

// Real space and self energy functions
void computeRealSpacePGPImpl(model::MCState& state, bool movement_only, bool store_in_residues = true);

double computeSelfEnergyPGPImpl(model::MCState& state, bool movement_only);

// Complete system energy evaluation functions
void computeSystemEnergyPGPImpl(model::MCState& state);

void computeMovementEnergyPGPImpl(model::MCState& state);

// <agent-hook:pgp_evaluator>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 