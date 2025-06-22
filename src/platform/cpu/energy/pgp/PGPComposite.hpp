#pragma once

#include "model/ModelModule.hpp"
#include "PGPCore.hpp"
#include "PGPGrid.hpp"
#include "PGPInterpolation.hpp"
#include "PGPPrecompute.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "PGPSystem.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Composite interface for PGP calculations
 * 
 * This module provides a unified interface for PGP energy calculations,
 * orchestrating all the individual PGP components (grid operations, interpolation,
 * precomputation, evaluation) to provide high-level energy calculation functions.
 */
class PGPComposite {
public:
    /**
     * @brief Initialize PGP with automatic parameter optimization
     * 
     * @param cutoff Real space cutoff distance
     * @param box Simulation box dimensions
     * @param alpha Ewald parameter
     * @param meshSize Grid dimensions
     * @param potentialCutoff Cutoff for potential calculation
     * @param potentialGridSize Precomputed potential grid dimensions
     * @param splineOrder B-spline order
     * @param tolerance Error tolerance
     */
    static void initialize(double cutoff, 
                         const double box[3],
                         double alpha,
                         const int meshSize[3],
                         double potentialCutoff,
                         const int potentialGridSize[3],
                         int splineOrder = 4,
                         double tolerance = 1e-5);

    /**
     * @brief Compute total system energy using PGP
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
     * @brief Validate PGP setup and parameters
     * 
     * @param state MC state
     * @return true if PGP is properly configured
     */
    static bool validateSetup(const model::MCState& state);

    /**
     * @brief Get PGP energy components breakdown
     * 
     * @param state MC state
     * @param realSpace Output: real space energy
     * @param reciprocal Output: reciprocal space energy (from precomputed grid)
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
     * @brief Check if PGP is properly initialized
     * 
     * @return true if PGP parameters are initialized and ready
     */
    static bool isInitialized();

    /**
     * @brief Reset PGP state for new calculation
     */
    static void reset();

    /**
     * @brief Enable or disable debug output
     * 
     * @param enable Whether to enable debug logging
     */
    static void setDebugMode(bool enable);

    /**
     * @brief Precompute potential grids
     * 
     * @param state MC state
     * @param fixedOnly Whether to only consider fixed residues
     */
    static void precomputeGrids(model::MCState& state, bool fixedOnly = true);

private:
    // Private helper functions for internal coordination
    static void validateSystemProperties(const model::MCState& state);
    static void checkGridCompatibility(const double box[3]);
    static void initializeGrids();
};

// Convenience functions that delegate to PGPComposite
// These maintain backward compatibility with existing code

/**
 * @brief Initialize PGP parameters (convenience function)
 */
inline void initializePGPParameters(double cutoff, const double box[3],
                                  double alpha,
                                  const int meshSize[3], 
                                  double potentialCutoff,
                                  const int potentialGridSize[3],
                                  int splineOrder = 4,
                                  double tolerance = 1e-5) {
    PGPComposite::initialize(cutoff, box, alpha, meshSize, 
                           potentialCutoff, potentialGridSize,
                           splineOrder, tolerance);
}

/**
 * @brief Compute system energy using PGP (convenience function)
 */
inline void computeSystemEnergyPGP(model::MCState& state) {
    PGPComposite::computeSystemEnergy(state);
}

/**
 * @brief Compute movement energy using PGP (convenience function)
 */
inline void computeMovementEnergyPGP(model::MCState& state) {
    PGPComposite::computeMovementEnergy(state);
}

/**
 * @brief Precompute grid potential using PGP (convenience function)
 */
inline void precomputeGridPotential(model::MCState& state, bool fixedOnly) {
    PGPComposite::precomputeGrids(state, fixedOnly);
}

// <agent-hook:pgp_composite>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 