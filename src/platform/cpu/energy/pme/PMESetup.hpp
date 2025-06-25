#pragma once

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Function Reference Guide:
 * - Parameter setting: setPMEParameters, autoAdjustPMEParameters
 * - Initialization: initializePMEParameters, initializePMETables, initializePMEBsplines
 * - Box management: setPMEBox
 * - Error estimation: estimate*Error functions
 */

// Constants for PME calculation
static const int DEFAULT_SPLINE_ORDER = 4;  // B-spline order for PME

// === Parameter Setting Functions ===

/**
 * @brief Set PME parameters
 * 
 * @param alpha Ewald separation parameter
 * @param meshSize Grid dimensions
 * @param splineOrder B-spline order (typically 4-6)
 * @param tolerance Error tolerance
 */
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder = DEFAULT_SPLINE_ORDER, double tolerance = 1e-5);

/**
 * @brief Auto-adjust PME parameters based on system properties
 * 
 * @param error_tolerance Target error tolerance
 * @param cutoff_distance Real-space cutoff distance
 * @param box Box dimensions
 */
void autoAdjustPMEParameters(double error_tolerance, double cutoff_distance, const double box[3]);

/**
 * @brief Initialize PME parameters with automatic optimization
 * 
 * @param cutoff Real space cutoff distance
 * @param box Simulation box dimensions
 * @param alpha Ewald separation parameter (auto-calculated if <= 0)
 * @param meshSize Grid dimensions (auto-calculated if null)
 * @param splineOrder B-spline order (typically 4-6)
 * @param tolerance Error tolerance
 */
void initializePMEParameters(double cutoff, const double box[3], 
                           double alpha, 
                           const int* meshSize,
                           int splineOrder,
                           double tolerance);

// === Table and Grid Initialization Functions ===

/**
 * @brief Initialize lookup tables for erfc and scaling functions
 * 
 * @param cutoff Cutoff distance
 */
void initializePMETables(double cutoff);

/**
 * @brief Set box dimensions for PME calculations
 * 
 * @param newBox Box dimensions [3]
 */
void setPMEBox(const double newBox[3]);

/**
 * @brief Initialize B-splines for PME
 */
void initializePMEBsplines();

// === Approximation Functions ===

/**
 * @brief Approximate erfc function using lookup table
 * 
 * @param r Distance
 * @return Approximate erfc value
 */
double erfcApproximate(double r);

/**
 * @brief Approximate electrostatic scaling function using lookup table
 * 
 * @param r Distance  
 * @return Scaling factor
 */
double ewaldScaleApproximate(double r);

// === Error Estimation Functions ===

/**
 * @brief Estimate real space error
 * 
 * @return Real space error estimate
 */
double estimatePMERealSpaceError();

/**
 * @brief Estimate reciprocal space error
 * 
 * @param box Box dimensions
 * @return Reciprocal space error estimate
 */
double estimatePMEReciprocalSpaceError(const double box[3]);

/**
 * @brief Estimate total PME error
 * 
 * @param box Box dimensions
 * @return Total error estimate
 */
double estimatePMETotalError(const double box[3]);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 