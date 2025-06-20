#pragma once

#include "model/montecarlo.hpp"
#include "PMECore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Reciprocal space calculations for PME
 * 
 * This module handles the reciprocal space (Fourier space) energy calculations
 * in the PME algorithm. This includes FFT operations and energy computation
 * from the transformed grid.
 */

/**
 * @brief Perform forward FFT on the PME grid
 * 
 * Transforms the real space charge distribution to reciprocal space
 * using the custom FFT implementation.
 */
void performFFTForward();

/**
 * @brief Perform backward FFT on the PME grid
 * 
 * Transforms reciprocal space data back to real space.
 * Used for force calculations (not needed for energy-only MC).
 */
void performFFTBackward();

/**
 * @brief Compute energy from the PME grid after FFT
 * 
 * This function calculates the reciprocal space contribution to the
 * electrostatic energy from the FFT-transformed grid.
 * 
 * @param energy Output variable for computed energy
 * @param box Simulation box dimensions
 */
void computeEnergyFromGrid(double& energy, const double box[3]);

/**
 * @brief Compute total reciprocal space energy using PME
 * 
 * High-level function that orchestrates the reciprocal space calculation:
 * charge spreading, FFT, and energy computation.
 * 
 * @param state MC state containing system information
 * @return Reciprocal space energy contribution
 */
double computeReciprocalEnergy(model::MCState& state);

/**
 * @brief Apply structure factor corrections in reciprocal space
 * 
 * @param grid PME grid in reciprocal space
 * @param params PME parameters
 * @param box Box dimensions
 */
void applyStructureFactorCorrections(std::vector<std::complex<double>>& grid,
                                   const PMEParams& params,
                                   const double box[3]);

/**
 * @brief Calculate reciprocal lattice vectors
 * 
 * @param box Real space box vectors
 * @param recipBox Output reciprocal space vectors
 */
void calculateReciprocalVectors(const double box[3], double recipBox[3][3]);

/**
 * @brief Validate reciprocal space calculation parameters
 * 
 * @param params PME parameters to validate
 * @param box Box dimensions
 * @return true if parameters are valid for reciprocal space calculation
 */
bool validateReciprocalParameters(const PMEParams& params, const double box[3]);

// <agent-hook:pme_reciprocal>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 