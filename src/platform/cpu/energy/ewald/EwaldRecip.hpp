#pragma once

#include "EwaldCore.hpp"
#include "model/ModelModule.hpp"
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Ewald reciprocal space calculations module
 *
 * This module handles the reciprocal space part of Ewald summation, including:
 * - Structure factor calculations
 * - K-space energy summation
 * - Fourier transform operations
 * - Convergence optimization
 */

/**
 * @brief Calculate reciprocal space energy using Ewald summation
 *
 * Computes the reciprocal space contribution to the total electrostatic energy
 * using the formula: E_recip = 4π/V * Σ_k (exp(-k²/4α²)/k²) * |S(k)|²
 *
 * @param state MC state containing system information
 * @param movement_only Whether to calculate only for moving residues
 * @return Reciprocal space energy contribution
 */
double computeReciprocalEnergy(model::MCState& state, bool movement_only);

/**
 * @brief Calculate structure factor for a given k-vector
 *
 * Computes S(k) = Σ_j q_j * exp(ik·r_j) for the given k-vector
 *
 * @param state MC state
 * @param kx K-vector x component
 * @param ky K-vector y component
 * @param kz K-vector z component
 * @param movement_only Whether to include only moving atoms
 * @return Complex structure factor
 */
std::complex<double> calculateStructureFactor(const model::MCState& state,
                                             double kx, double ky, double kz,
                                             bool movement_only);

/**
 * @brief Optimize k-vector summation limits
 *
 * Determines optimal kmax values for each dimension based on
 * convergence criteria and computational efficiency.
 *
 * @param box Simulation box dimensions
 * @param alpha Ewald separation parameter
 * @param tolerance Desired accuracy
 * @param kmax Output array for optimized kmax values
 */
void optimizeKVectorLimits(const double box[3],
                          double alpha,
                          double tolerance,
                          int kmax[3]);

/**
 * @brief Validate reciprocal space calculation parameters
 *
 * @param state MC state
 * @return true if parameters are valid for reciprocal space calculation
 */
bool validateReciprocalSpaceParameters(const model::MCState& state);

/**
 * @brief Get reciprocal space convergence information
 *
 * @param state MC state
 * @param numKVectors Output: total number of k-vectors used
 * @param maxKVector Output: maximum k-vector magnitude
 * @param estimatedError Output: estimated convergence error
 */
void getReciprocalSpaceInfo(const model::MCState& state,
                          int& numKVectors,
                          double& maxKVector,
                          double& estimatedError);

// <agent-hook:ewald_reciprocal>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
