#include "energyPGP.hpp"

/**
 * @file energyPGP.cpp
 * @brief Implementation of Precomputed Grid-Potential Particle Mesh Ewald for Monte Carlo (PGP-PME-MC)
 * 
 * This file now serves as a facade that includes all PGP submodules.
 * The actual implementations have been split into separate files for better maintainability:
 * - PGPGrid.cpp: Grid initialization and management
 * - PGPPrecompute.cpp: Parameter setup and potential precomputation  
 * - PGPInterpolation.cpp: Energy interpolation calculations
 * - PGPEvaluator.cpp: Real space, self energy, and system energy evaluation
 */

namespace pygcmc {
namespace platform {
namespace cpu {

// Initialize global PGP parameters
PGPParams pgp_params;

// All implementations are now in the submodules
// <agent-hook:pgp_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 