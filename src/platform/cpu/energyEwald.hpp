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
static const int NUM_TABLE_POINTS = 2048;
static const float TWO_OVER_SQRT_PI = 2.0f/std::sqrt(M_PI);

// Ewald parameters
struct EwaldParams {
    float alpha{1.0f};     // Ewald separation parameter (nm^-1)
    int kmax[3]{6,6,6};    // Maximum reciprocal space wave vectors
    float tolerance{1e-5f}; // Precision control
    bool initialized{false};
    float cutoff{0.0f};    // Real space cutoff distance
    
    // Tables for optimized calculation
    std::vector<float> erfcTable;      // Table for erfc values
    std::vector<float> ewaldScaleTable; // Table for complete Ewald scale factor
    float ewaldDX;                     // Table spacing
    float ewaldDXInv;                  // Inverse of table spacing
    float erfcDXInv;                   // Inverse of table spacing for erfc
    
    // Exp(ikr) tables for reciprocal space optimization
    std::vector<std::complex<float>> expIkrTable;  // Table for exp(ikr)
    std::vector<std::complex<float>> expIkrXY;     // Temporary storage for xy plane
    int maxK;                          // Maximum k value for tables
    
    // Methods for table management
    void initializeTables(float cutoff);
    void initializeExpIkrTable(int numAtoms);
    float erfcApprox(float r) const;
    float ewaldScaleApprox(float r) const;
    
    // Error estimation methods
    float estimateRealSpaceError() const {
        return std::erfc(alpha * cutoff);
    }
    
    float estimateReciprocalSpaceError(const float box[3]) const {
        float minBoxSize = std::min(box[0], std::min(box[1], box[2]));
        float minKmax = std::min(kmax[0], std::min(kmax[1], kmax[2]));
        float error = minKmax * std::sqrt(alpha * minBoxSize) / 20.0f;
        error *= std::exp(-(M_PI * minKmax / (alpha * minBoxSize)) * 
                         (M_PI * minKmax / (alpha * minBoxSize)));
        return error;
    }
    
    float estimateTotalError(const float box[3]) const {
        return estimateRealSpaceError() + estimateReciprocalSpaceError(box);
    }
};

// Global Ewald parameters
extern EwaldParams ewald_params;

// Function to set Ewald parameters
void setEwaldParameters(float alpha, const int kmax[3], float tolerance = 1e-5f);

// Automatic parameter selection
void autoAdjustParameters(float error_tolerance, float cutoff_distance, const float box[3]);

// Ewald method interfaces
void computeSystemEnergyEwald(model::MCState& state);
void computeMovementEnergyEwald(model::MCState& state);

// Internal calculation methods
float computeReciprocalEnergy(model::MCState& state, bool movement_only);
float computeSelfEnergy(model::MCState& state, bool movement_only);
std::pair<float, float> calcPairEnergyEwald(float r2, float sigma, float eps, float q1, float q2);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 