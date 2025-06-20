#pragma once

#include "model/montecarlo.hpp"
#include "EwaldCore.hpp"
#include "EwaldReal.hpp"
#include "EwaldRecip.hpp"
#include "EwaldSelf.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Composite interface for Ewald calculations
 * 
 * This module provides a unified interface for Ewald energy calculations,
 * orchestrating all the individual Ewald components (real space, reciprocal space,
 * self-energy) to provide high-level energy calculation functions.
 */
class EwaldComposite {
public:
    /**
     * @brief Initialize Ewald with automatic parameter optimization
     * 
     * @param cutoff Real space cutoff distance
     * @param box Simulation box dimensions
     * @param alpha Ewald parameter (auto-calculated if <= 0)
     * @param tolerance Error tolerance
     */
    static void initialize(double cutoff, 
                         const double box[3], 
                         double alpha = 0.0,
                         double tolerance = 1e-5);

    /**
     * @brief Compute total system energy using Ewald summation
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
     * @brief Validate Ewald setup and parameters
     * 
     * @param state MC state
     * @return true if Ewald is properly configured
     */
    static bool validateSetup(const model::MCState& state);

    /**
     * @brief Get Ewald energy components breakdown
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
     * @brief Check if Ewald is properly initialized
     * 
     * @return true if Ewald parameters are initialized and ready
     */
    static bool isInitialized();

    /**
     * @brief Reset Ewald state for new calculation
     */
    static void reset();

    /**
     * @brief Get convergence information for all components
     * 
     * @param state MC state
     * @param realSpaceError Output: estimated real space error
     * @param reciprocalError Output: estimated reciprocal space error
     * @param totalError Output: estimated total error
     */
    static void getConvergenceInfo(const model::MCState& state,
                                 double& realSpaceError,
                                 double& reciprocalError,
                                 double& totalError);

private:
    // Private helper functions for internal coordination
    static void validateSystemProperties(const model::MCState& state);
    static void checkParameterConsistency();
};

// Convenience functions that delegate to EwaldComposite
// These maintain backward compatibility with existing code

/**
 * @brief Initialize Ewald parameters (convenience function)
 */
inline void initializeEwaldParameters(double cutoff, const double box[3], 
                                    double alpha = 0.0, double tolerance = 1e-5) {
    EwaldComposite::initialize(cutoff, box, alpha, tolerance);
}

/**
 * @brief Compute system energy using Ewald (convenience function)
 */
inline void computeSystemEnergyEwald(model::MCState& state) {
    EwaldComposite::computeSystemEnergy(state);
}

/**
 * @brief Compute movement energy using Ewald (convenience function)
 */
inline void computeMovementEnergyEwald(model::MCState& state) {
    EwaldComposite::computeMovementEnergy(state);
}

/**
 * @brief Calculate pair energy using Ewald (convenience function)
 */
inline std::pair<double, double> calcPairEnergyEwald(double r2, double sigma, double eps, 
                                                   double q1, double q2, 
                                                   const model::MCInfo& info, 
                                                   bool is_excluded = false) {
    return calcPairEnergyEwaldRealSpace(r2, sigma, eps, q1, q2, info, is_excluded);
}

// <agent-hook:ewald_composite>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 