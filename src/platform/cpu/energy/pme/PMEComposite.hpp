#pragma once

#include "model/montecarlo.hpp"
#include "PMECore.hpp"
#include "PMESetup.hpp"
#include "PMEInterface.hpp"
#include "PMEFFT.hpp"
#include "PMESpline.hpp"
#include "PMEGrid.hpp"
#include "PMERecip.hpp"
#include "PMEReal.hpp"
#include "PMESelf.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Composite interface for PME calculations
 * 
 * This module provides a unified interface for PME energy calculations,
 * orchestrating all the individual PME components (real space, reciprocal space,
 * self-energy, etc.) to provide high-level energy calculation functions.
 */
class PMEComposite {
public:
    /**
     * @brief Initialize PME with automatic parameter optimization
     * 
     * @param cutoff Real space cutoff distance
     * @param box Simulation box dimensions
     * @param tolerance Error tolerance
     * @param alpha Ewald parameter (auto-calculated if <= 0)
     * @param meshSize Grid dimensions (auto-calculated if null)
     * @param splineOrder B-spline order
     */
    static void initialize(double cutoff, 
                         const double box[3], 
                         double tolerance = 1e-5,
                         double alpha = 0.0,
                         const int* meshSize = nullptr,
                         int splineOrder = 4);

    /**
     * @brief Compute total system energy using PME
     * 
     * Calculates all components: real space + reciprocal space + self energy + VdW
     * 
     * @param state MC state containing system information
     */
    static void computeSystemEnergy(model::MCState& state);

    /**
     * @brief Compute energy for moving residues only
     * 
     * Optimized calculation for Monte Carlo moves that only affect
     * a subset of the system.
     * 
     * @param state MC state
     */
    static void computeMovementEnergy(model::MCState& state);

    /**
     * @brief Validate PME setup and parameters
     * 
     * @param state MC state
     * @return true if PME is properly configured
     */
    static bool validateSetup(const model::MCState& state);

    /**
     * @brief Get PME energy components breakdown
     * 
     * @param state MC state
     * @param realSpace Output: real space energy
     * @param reciprocal Output: reciprocal space energy  
     * @param selfEnergy Output: self energy
     * @param vdw Output: van der Waals energy
     * @param total Output: total energy
     */
    static void getEnergyBreakdown(const model::MCState& state,
                                 double& realSpace,
                                 double& reciprocal,
                                 double& selfEnergy,
                                 double& vdw,
                                 double& total);

    /**
     * @brief Check if PME is properly initialized
     * 
     * @return true if PME parameters are initialized and ready
     */
    static bool isInitialized();

    /**
     * @brief Reset PME state for new calculation
     */
    static void reset();

    /**
     * @brief Enable or disable debug output
     * 
     * @param enable Whether to enable debug logging
     */
    static void setDebugMode(bool enable);

private:
    // Private helper functions for internal coordination
    static void validateSystemProperties(const model::MCState& state);
    static void checkGridCompatibility(const double box[3]);
};

// Convenience functions that delegate to PMEComposite
// These maintain backward compatibility with existing code

/**
 * @brief Initialize PME parameters (convenience function)
 * 
 * @param cutoff Cutoff distance
 * @param box Box dimensions
 * @param alpha Ewald separation parameter (if <=0, automatically calculated)
 * @param meshSize Mesh dimensions for PME (if NULL, automatically calculated)
 * @param splineOrder B-spline order (typically 4-6)
 * @param tolerance Precision control parameter
 */
inline void initializePMEParameters(double cutoff, const double box[3], 
                                  double alpha = 0.0, 
                                  const int* meshSize = nullptr,
                                  int splineOrder = DEFAULT_SPLINE_ORDER,
                                  double tolerance = 1e-5) {
    // Set box dimensions - ensure B-spline initialization uses correct volume
    pme_params.setBox(box);
    
    // If alpha is not specified, calculate the optimal value
    if (alpha <= 0.0) {
        // Auto-adjust parameters includes calling initializeTables and initializeBsplines
        autoAdjustPMEParameters(tolerance, cutoff, box);
    } else {
        // Use the specified parameters
        int mSize[3] = {64, 64, 64}; // Default value
        
        // If mesh size is provided, use it
        if (meshSize != nullptr) {
            mSize[0] = meshSize[0];
            mSize[1] = meshSize[1];
            mSize[2] = meshSize[2];
        }
        
        // Set parameters and initialize tables
        setPMEParameters(alpha, mSize, splineOrder, tolerance);
        pme_params.initializeTables(cutoff);
        pme_params.initializeBsplines();
    }
    
    // Log the final parameters for debugging
    platform::log(LogLevel::INFO, "PME parameters initialized: alpha = ", pme_params.alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder,
                 ", initialized = ", pme_params.initialized);
}

/**
 * @brief Compute system energy using PME (convenience function)
 */
inline void computeSystemEnergyPME(model::MCState& state) {
    PMEComposite::computeSystemEnergy(state);
}

/**
 * @brief Compute movement energy using PME (convenience function)
 */
inline void computeMovementEnergyPME(model::MCState& state) {
    PMEComposite::computeMovementEnergy(state);
}

// <agent-hook:pme_composite>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 