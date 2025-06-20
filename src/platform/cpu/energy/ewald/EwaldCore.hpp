#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "platform/cpu/energyCommon.hpp"
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Ewald summation algorithm parameter structure
 * 
 * This structure contains all parameters and data structures required by the Ewald algorithm.
 * Ewald summation is a method for computing long-range electrostatic interactions in 
 * periodic systems by splitting into real space and reciprocal space contributions.
 */
struct EwaldParams {
    double alpha{1.0};                 // Ewald separation parameter
    int kmax[3]{15,15,15};             // Maximum reciprocal space wave vectors
    double tolerance{1e-5f};           // Error tolerance
    bool initialized{false};           // Whether parameters have been initialized
    double cutoff{0.0};               // Real space cutoff distance
    
    // Lookup tables for performance optimization
    std::vector<double> erfcTable;         // erfc function lookup table
    std::vector<double> ewaldScaleTable;   // Ewald scaling factor lookup table
    double ewaldDX;                        // Table step size
    double ewaldDXInv;                     // Inverse of table step size
    double erfcDXInv;                      // Inverse of erfc table step size
    
    // Reciprocal space optimization with exp(ikr) tables
    std::vector<std::complex<double>> expIkrTable;  // Pre-computed exp(ikr) values
    std::vector<std::complex<double>> expIkrXY;     // XY plane exp(ikr) values
    int maxK;                              // Maximum k value for pre-computation
    
    // Core parameter management methods
    void initializeTables(double cutoff);
    void initializeExpIkrTable(int numAtoms);
    double erfcApprox(double r) const;
    double ewaldScaleApprox(double r) const;
    
    // Error estimation methods
    double estimateRealSpaceError() const;
    double estimateReciprocalSpaceError(const double box[3]) const;
    double estimateTotalError(const double box[3]) const;
};

// Global Ewald parameters instance
extern EwaldParams ewald_params;

// Mathematical constants
extern const double TWO_PI;
extern const double SQRT_PI;

// Core initialization and parameter setting functions
void setEwaldParameters(double alpha, const int kmax[3], double tolerance = 1e-5);
void autoAdjustParameters(double error_tolerance, double cutoff_distance, const double box[3]);

/**
 * @brief Initialize Ewald parameters with automatic optimization
 * 
 * @param cutoff Real space cutoff distance
 * @param box Simulation box dimensions
 * @param alpha Ewald separation parameter (auto-calculated if <= 0)
 * @param tolerance Error tolerance
 */
void initializeEwaldParameters(double cutoff, const double box[3], 
                             double alpha, double tolerance);

// High-level energy calculation interfaces
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);

// Component energy calculation functions
double computeReciprocalEnergy(model::MCState& state, bool movement_only);
double computeSelfEnergy(model::MCState& state, bool movement_only);
void computeRealSpaceEwald(model::MCState& state, bool movement_only, bool store_in_residues);

// Pair energy calculation for Ewald
std::pair<double, double> calcPairEnergyEwald(double r2, double sigma, double eps, 
                                            double q1, double q2, 
                                            const model::MCInfo& info, 
                                            bool is_excluded);

// <agent-hook:ewald_core>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 