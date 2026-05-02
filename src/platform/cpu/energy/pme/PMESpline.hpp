#pragma once

#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief B-spline interpolation functions for PME
 *
 * This module provides functions for computing B-spline coefficients and moduli
 * used in the PME charge interpolation process. B-splines provide smooth
 * interpolation between grid points for accurate charge distribution.
 */

/**
 * @brief Compute B-spline coefficients for a given fractional position
 *
 * @param fractional Fractional position (0-1)
 * @param order Spline order (typically 4)
 * @param coefficients Output vector for coefficients
 */
void computeBSplineCoefficients(double fractional, int order, std::vector<double>& coefficients);

/**
 * @brief Initialize B-spline moduli for PME calculations
 *
 * This function computes the B-spline moduli needed for the reciprocal space
 * energy calculation. The moduli are computed for each grid dimension.
 *
 * @param params PME parameters structure containing grid dimensions and spline order
 */
void initializeBSplineModuli(struct PMEParams& params);

/**
 * @brief Compute B-spline interpolation weights for charge spreading
 *
 * @param position Particle position in fractional coordinates
 * @param order Spline order
 * @param weights Output array for interpolation weights
 */
void computeInterpolationWeights(const double position[3], int order, double* weights);

/**
 * @brief Validate B-spline coefficients for consistency
 *
 * @param coefficients B-spline coefficients to validate
 * @param order Spline order
 * @return true if coefficients are valid, false otherwise
 */
bool validateBSplineCoefficients(const std::vector<double>& coefficients, int order);

// <agent-hook:pme_spline>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
