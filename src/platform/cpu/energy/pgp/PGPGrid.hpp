#pragma once

#include "PGPCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Grid management functions for PGP-PME
 *
 * This module handles the initialization and management of the potential grid
 * used in the Precomputed Grid-Potential Particle Mesh Ewald algorithm.
 */

// Grid initialization and management functions
void initializePotentialGridImpl();

// <agent-hook:pgp_grid>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
