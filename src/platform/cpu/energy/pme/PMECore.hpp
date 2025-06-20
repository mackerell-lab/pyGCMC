#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "platform/cpu/energyCommon.hpp"
#include <array>
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for PME calculation
static const int DEFAULT_SPLINE_ORDER = 4;  // B-spline order for PME

/**
 * @brief PME (Particle Mesh Ewald) algorithm parameter structure
 * 
 * This structure contains all parameters and data structures required by the PME algorithm.
 * PME is a widely used method for computing long-range electrostatic interactions in 
 * periodic systems by decomposing the calculation into real space and reciprocal space parts.
 */
struct PMEParams {
    double alpha{1.0};                 // Ewald separation parameter
    int meshSize[3]{64,64,64};         // Default mesh size for PME
    double tolerance{1e-5f};           // Error tolerance for PME accuracy
    bool initialized{false};           // Whether parameters have been initialized
    double cutoff{0.0};               // Real space cutoff distance
    double epsilon_r{1.0};            // Relative dielectric constant
    double box[3]{1.0, 1.0, 1.0};     // Box dimensions, default to unit box
    
    // Lookup tables for performance optimization
    std::vector<double> erfcTable;         // erfc function lookup table
    std::vector<double> ewaldScaleTable;   // Ewald scaling factor lookup table
    double ewaldDX;                        // Table step size
    double ewaldDXInv;                     // Inverse of table step size
    double erfcDXInv;                      // Inverse of erfc table step size

    // B-spline parameters
    int splineOrder{DEFAULT_SPLINE_ORDER}; // B-spline order (typically 4-6)
    std::vector<double> bsplineModuli[3];  // Bspline moduli for the 3 dimensions
    
    // Grid data structures
    std::vector<std::complex<double>> pmeGrid;   // Main PME grid for FFT operations
    std::vector<double> pmeCharge;              // Charge grid (if needed)
    
    // Core parameter management methods
    void initializeTables(double cutoff);
    void initializeBsplines();
    void setBox(const double newBox[3]);
    double erfcApprox(double r) const;
    double ewaldScaleApprox(double r) const;
    
    // Error estimation methods
    double estimateRealSpaceError() const;
    double estimateReciprocalSpaceError(const double box[3]) const;
    double estimateTotalError(const double box[3]) const;
};

// Global PME parameters instance
extern PMEParams pme_params;

// Core initialization and parameter setting functions
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder = DEFAULT_SPLINE_ORDER, double tolerance = 1e-5);
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

// High-level energy calculation interfaces
void computeSystemEnergyPME(model::MCState& state);
void computeMovementEnergyPME(model::MCState& state);

// Component energy calculation functions
double computeReciprocalPME(model::MCState& state);
double computeSelfEnergyPME(model::MCState& state, bool movement_only);
void computeRealSpacePME(model::MCState& state, bool movement_only, bool store_in_residues = true);

// Pair energy calculation for PME
std::pair<double, double> calcPairEnergyPME(double r2, double sigma, double eps, 
                                          double q1, double q2, 
                                          const model::MCInfo& info, 
                                          bool is_excluded = false);

// <agent-hook:pme_core>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 