#ifndef EWALDMOVEMENTENERGY_HPP
#define EWALDMOVEMENTENERGY_HPP

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute energy for moving residues only
 * 
 * @param state MC state
 */
void computeMovementEwaldEnergy(model::MCState& state);

/**
 * @brief Initialize Ewald with automatic parameter optimization
 * 
 * @param cutoff Cutoff distance
 * @param box Box dimensions
 * @param alpha Alpha parameter (0 for auto)
 * @param tolerance Error tolerance
 */
void initializeEwald(double cutoff, const double box[3], double alpha = 0.0, double tolerance = 1e-5);

/**
 * @brief Check if Ewald is properly initialized
 * 
 * @return true if initialized
 */
bool isEwaldInitialized();

/**
 * @brief Reset Ewald state for new calculation
 */
void resetEwald();

/**
 * @brief Get convergence information for all components
 * 
 * @param state MC state
 * @param realSpaceError Real space error estimate
 * @param reciprocalError Reciprocal space error estimate
 * @param totalError Total error estimate
 */
void getEwaldConvergenceInfo(const model::MCState& state,
                            double& realSpaceError,
                            double& reciprocalError,
                            double& totalError);

/**
 * @brief Check parameter consistency
 */
void checkEwaldParameterConsistency();

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // EWALDMOVEMENTENERGY_HPP 