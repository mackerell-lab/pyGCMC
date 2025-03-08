// src/platform/cpu/energyEwald.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for Ewald calculation already defined in energyCommon.hpp
// static const int NUM_TABLE_POINTS = 20000;
// static const double TWO_OVER_SQRT_PI = 2.0/std::sqrt(M_PI);

// Ewald parameters structure
struct EwaldParams {
    double alpha{1.0};     // Changed to double for higher precision
    int kmax[3]{15,15,15}; // Increased from 6,6,6 to 15,15,15 for better convergence
    double tolerance{1e-5f};
    bool initialized{false};
    double cutoff{0.0};    // Changed to double
    
    // Lookup tables for optimization
    std::vector<double> erfcTable;      // Changed to double
    std::vector<double> ewaldScaleTable;
    double ewaldDX;                     // Changed to double
    double ewaldDXInv;
    double erfcDXInv;
    
    // Reciprocal space optimization with exp(ikr) tables
    std::vector<std::complex<double>> expIkrTable;  // Changed to double
    std::vector<std::complex<double>> expIkrXY;
    int maxK;
    
    // Table management methods
    void initializeTables(double cutoff);
    void initializeExpIkrTable(int numAtoms);
    double erfcApprox(double r) const;
    double ewaldScaleApprox(double r) const;
    
    // Error estimation methods
    double estimateRealSpaceError() const {
        return std::erfc(alpha * cutoff);
    }
    
    double estimateReciprocalSpaceError(const double box[3]) const {
        double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
        double minKmax = std::min(kmax[0], std::min(kmax[1], kmax[2]));
        double error = minKmax * std::sqrt(alpha * minBoxSize) / 20.0;
        error *= std::exp(-(M_PI * minKmax / (alpha * minBoxSize)) * 
                         (M_PI * minKmax / (alpha * minBoxSize)));
        return error;
    }
    
    double estimateTotalError(const double box[3]) const {
        return estimateRealSpaceError() + estimateReciprocalSpaceError(box);
    }
};

// Global Ewald parameters
extern EwaldParams ewald_params;

// Function declarations
void setEwaldParameters(double alpha, const int kmax[3], double tolerance = 1e-5);
void autoAdjustParameters(double error_tolerance, double cutoff_distance, const double box[3]);

/**
 * @brief Initialize Ewald parameters
 * 
 * @param cutoff Cutoff distance
 * @param box Box dimensions
 * @param alpha Ewald separation parameter (if <=0, automatically calculated)
 * @param tolerance Precision control parameter
 */
inline void initializeEwaldParameters(double cutoff, const double box[3], 
                                     double alpha = 0.0, double tolerance = 1e-5) {
    // If alpha is not specified, calculate the optimal value
    if (alpha <= 0.0) {
        autoAdjustParameters(tolerance, cutoff, box);
    } else {
        // Use the specified alpha value
        int kmax[3] = {15, 15, 15}; // Default value
        setEwaldParameters(alpha, kmax, tolerance);
        ewald_params.initializeTables(cutoff);
    }
}

/**
 * @brief Use Ewald method to calculate system energy
 * 
 * @param state MC state
 */
void computeSystemEnergyEwald(model::MCState& state);

/**
 * @brief Use Ewald method to calculate energy of moving residue
 * 
 * @param state MC state
 */
void computeMovementEnergyEwald(model::MCState& state);

// Internal calculation function declarations
std::pair<double, double> calcPairEnergyEwald(double r2, double sigma, double eps, double q1, double q2, const model::MCInfo& info, bool is_excluded = false);
double computeReciprocalEnergy(model::MCState& state, bool movement_only);
double computeSelfEnergy(model::MCState& state, bool movement_only);
void computeRealSpaceEwald(model::MCState& state, bool movement_only, bool store_in_residues = true);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 