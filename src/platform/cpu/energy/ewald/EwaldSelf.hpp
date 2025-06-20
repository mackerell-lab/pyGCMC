#pragma once

#include "EwaldCore.hpp"
#include "model/montecarlo.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Ewald self-energy corrections module
 * 
 * This module handles the self-energy correction term in Ewald summation.
 * The self-energy arises from the artificial interaction of each charge with
 * itself in the reciprocal space sum and must be subtracted to get the correct
 * total electrostatic energy.
 */

/**
 * @brief Calculate self-energy correction for Ewald summation
 * 
 * Computes the self-energy correction term using the formula:
 * E_self = -α/√π * Σ_i q_i²
 * 
 * @param state MC state containing system information
 * @param movement_only Whether to calculate only for moving residues
 * @return Self-energy correction (negative value)
 */
double computeSelfEnergy(model::MCState& state, bool movement_only);

/**
 * @brief Calculate self-energy for a specific set of charges
 * 
 * @param charges Vector of charges to compute self-energy for
 * @param alpha Ewald separation parameter
 * @return Self-energy correction
 */
double calculateSelfEnergyForCharges(const std::vector<double>& charges, double alpha);

/**
 * @brief Get self-energy breakdown by residue
 * 
 * @param state MC state
 * @param selfEnergies Output vector of per-residue self-energies
 */
void getSelfEnergyBreakdown(const model::MCState& state, 
                           std::vector<double>& selfEnergies);

/**
 * @brief Validate self-energy calculation parameters
 * 
 * @param state MC state
 * @return true if parameters are valid for self-energy calculation
 */
bool validateSelfEnergyParameters(const model::MCState& state);

/**
 * @brief Calculate self-energy contribution for a single charge
 * 
 * @param charge Charge value
 * @param alpha Ewald separation parameter
 * @return Self-energy contribution for this charge
 */
inline double calculateSingleChargeSelfEnergy(double charge, double alpha) {
    return -COULOMB * alpha / std::sqrt(M_PI) * charge * charge;
}

// <agent-hook:ewald_self>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 