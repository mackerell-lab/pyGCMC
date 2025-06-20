#pragma once

#include "EwaldCore.hpp"
#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Ewald real space calculations module
 * 
 * This module handles the real space part of Ewald summation, including:
 * - Pair energy calculations with erfc(αr)/r terms
 * - Real space energy accumulation
 * - Cutoff-based neighbor list processing
 */

/**
 * @brief Calculate pair energy for Ewald real-space part
 * 
 * For normal pairs: erfc(αr)/r
 * For excluded pairs: -erf(αr)/r to compensate for reciprocal space
 * 
 * @param r2 Squared distance between particles
 * @param sigma Lennard-Jones sigma parameter
 * @param eps Lennard-Jones epsilon parameter  
 * @param q1 Charge of first particle
 * @param q2 Charge of second particle
 * @param info Simulation info containing cutoffs and parameters
 * @param is_excluded Whether this is an excluded pair interaction
 * @return Pair containing VdW energy and electrostatic energy
 */
std::pair<double, double> calcPairEnergyEwaldRealSpace(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool is_excluded);

/**
 * @brief Calculate real-space part of Ewald summation
 * 
 * Computes the real-space energy contribution using erfc(αr)/r terms
 * for all particle pairs within the cutoff distance.
 * 
 * @param state MC state containing system information
 * @param movement_only Whether to calculate only for moving residues
 * @param store_in_residues Whether to store energy in residue objects
 */
void computeRealSpaceEwald(model::MCState& state, 
                          bool movement_only, 
                          bool store_in_residues);

/**
 * @brief Validate real space calculation parameters
 * 
 * @param state MC state
 * @return true if parameters are valid for real space calculation
 */
bool validateRealSpaceParameters(const model::MCState& state);

/**
 * @brief Get real space energy breakdown by residue
 * 
 * @param state MC state
 * @param energies Output vector of per-residue real space energies
 */
void getRealSpaceEnergyBreakdown(const model::MCState& state, 
                               std::vector<double>& energies);

// <agent-hook:ewald_real_space>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 