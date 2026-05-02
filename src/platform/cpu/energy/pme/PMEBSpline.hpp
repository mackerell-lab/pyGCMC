#ifndef PMEBSPLINE_HPP
#define PMEBSPLINE_HPP

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Initialize B-splines for PME
 */
void initializePMEBsplines();

/**
 * @brief Set PME parameters
 *
 * @param alpha Ewald alpha parameter
 * @param meshSize Grid dimensions
 * @param splineOrder B-spline order
 * @param tolerance Error tolerance
 */
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder, double tolerance);

/**
 * @brief Auto-adjust PME parameters based on system properties
 *
 * @param error_tolerance Target error tolerance
 * @param cutoff_distance Real-space cutoff distance
 * @param box Box dimensions
 */
void autoAdjustPMEParameters(double error_tolerance, double cutoff_distance, const double box[3]);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEBSPLINE_HPP
