// src/platform/cpu/energyEwald.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include "energyDirect.hpp"  // For direct calculation methods
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for Ewald calculation
static const int NUM_TABLE_POINTS = 20000;  // Increased from 2048 for better precision
static const double TWO_OVER_SQRT_PI = 2.0/std::sqrt(M_PI);

// Ewald parameters
struct EwaldParams {
    double alpha{1.0};     // Changed to double for better precision
    int kmax[3]{15,15,15}; // Increased from 6,6,6 for better convergence
    double tolerance{1e-5f};
    bool initialized{false};
    double cutoff{0.0};    // Changed to double
    
    // Tables for optimized calculation
    std::vector<double> erfcTable;      // Changed to double
    std::vector<double> ewaldScaleTable;
    double ewaldDX;                     // Changed to double
    double ewaldDXInv;
    double erfcDXInv;
    
    // Exp(ikr) tables for reciprocal space optimization
    std::vector<std::complex<double>> expIkrTable;  // Changed to double
    std::vector<std::complex<double>> expIkrXY;
    int maxK;
    
    // Methods for table management
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
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);
std::pair<double, double> calcPairEnergyEwald(double r2, double sigma, double eps, double q1, double q2);
double computeReciprocalEnergy(model::MCState& state, bool movement_only);
double computeSelfEnergy(model::MCState& state, bool movement_only);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 