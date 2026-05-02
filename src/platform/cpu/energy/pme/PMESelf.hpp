#pragma once

#include "model/ModelModule.hpp"
#include "PMECore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Self-energy calculations for PME
 *
 * This module handles the self-energy correction terms in the PME algorithm.
 * The self-energy corrects for the spurious self-interaction that arises from
 * the reciprocal space calculation.
 */

/**
 * @brief Compute self-energy correction for PME
 *
 * The self-energy term removes the spurious self-interaction energy that
 * each particle has with itself in the reciprocal space calculation.
 *
 * @param state MC state containing particle charges
 * @param movement_only If true, only calculate for moving particles
 * @return Self-energy correction (negative value)
 */
double computeSelfEnergyPME(model::MCState& state, bool movement_only = false);

/**
 * @brief Calculate self-energy for a single particle
 *
 * @param charge Particle charge
 * @param alpha PME alpha parameter
 * @return Self-energy contribution for this particle
 */
double calculateParticleSelfEnergy(double charge, double alpha);

/**
 * @brief Calculate total system charge for neutrality check
 *
 * @param state MC state
 * @param movement_only If true, only consider moving particles
 * @return Total system charge
 */
double calculateTotalSystemCharge(const model::MCState& state, bool movement_only = false);

/**
 * @brief Validate system neutrality for PME calculations
 *
 * PME requires a neutral system for proper convergence. This function
 * checks if the total system charge is sufficiently close to zero.
 *
 * @param totalCharge Total system charge
 * @param tolerance Tolerance for neutrality check
 * @return true if system is sufficiently neutral
 */
bool validateSystemNeutrality(double totalCharge, double tolerance = 1e-6);

/**
 * @brief Get self-energy prefactor
 *
 * @param alpha PME alpha parameter
 * @return Self-energy prefactor (-alpha/sqrt(pi))
 */
inline double getSelfEnergyPrefactor(double alpha) {
    return -alpha / std::sqrt(M_PI);
}

// <agent-hook:pme_self>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
