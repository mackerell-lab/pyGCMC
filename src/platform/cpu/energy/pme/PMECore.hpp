#pragma once

#include "model/ModelModule.hpp"
#include "platform/platform.hpp"
#include "../common/EnergyConstants.hpp"
#include "../common/EnergyUtils.hpp"
#include "PMESetup.hpp"
#include <array>
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Function Reference Guide:
 * - Parameter structure: PMEParams struct
 * - Setup and initialization: PMESetup.hpp
 * - High-level interfaces: PMEInterface.hpp
 * - Real space: PMERealSpace.hpp  
 * - Self energy: PMESelf.hpp
 * - Reciprocal space: PMEReciprocal.hpp
 */

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
    
    // Legacy member functions - kept for compatibility with existing code
    void initializeTables(double cutoff) {
        // Redirect to standalone function
        initializePMETables(cutoff);
    }
    
    void setBox(const double newBox[3]) {
        // Redirect to standalone function
        setPMEBox(newBox);
    }
    
    void initializeBsplines() {
        // Redirect to standalone function
        initializePMEBsplines();
    }
    
    double erfcApprox(double r) const {
        // Redirect to standalone function
        return erfcApproximate(r);
    }
    
    double ewaldScaleApprox(double r) const {
        // Redirect to standalone function
        return ewaldScaleApproximate(r);
    }
    
    double estimateRealSpaceError() const {
        // Redirect to standalone function
        return estimatePMERealSpaceError();
    }
    
    double estimateReciprocalSpaceError(const double box[3]) const {
        // Redirect to standalone function
        return estimatePMEReciprocalSpaceError(box);
    }
    
    double estimateTotalError(const double box[3]) const {
        // Redirect to standalone function
        return estimatePMETotalError(box);
    }
};

// Global PME parameters instance - now managed by smart pointer
// See PMEGlobal.hpp for access functions
// extern PMEParams pme_params; // DEPRECATED - use getPMEParams() instead

// Clear all PME state - used for testing to prevent cross-test contamination
void clearPMEState();

// <agent-hook:pme_core>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 