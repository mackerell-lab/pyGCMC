#pragma once

#include "model/ModelModule.hpp"
#include "PMECore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Real space calculations for PME
 * 
 * This module handles the real space (short-range) electrostatic interactions
 * in the PME algorithm. It computes the erfc-damped Coulomb interactions
 * within the cutoff distance.
 */

/**
 * @brief Compute real space PME energy for the entire system
 * 
 * @param state MC state containing particle positions and charges
 * @param movement_only If true, only calculate for moving residues
 * @param store_in_residues If true, store energy contributions in residue structures
 */
void computeRealSpaceEnergy(model::MCState& state, 
                          bool movement_only = false, 
                          bool store_in_residues = true);

/**
 * @brief Calculate pair energy using PME real space formula
 * 
 * Computes the short-range electrostatic and LJ interactions between
 * two particles using the PME real space formulation.
 * 
 * @param r2 Squared distance between particles
 * @param sigma LJ sigma parameter
 * @param eps LJ epsilon parameter
 * @param q1, q2 Particle charges
 * @param info MC info containing simulation parameters
 * @param is_excluded Whether this is an excluded pair interaction
 * @return Pair of (LJ energy, electrostatic energy)
 */
std::pair<double, double> calculatePairEnergyPME(double r2, 
                                                double sigma, 
                                                double eps, 
                                                double q1, 
                                                double q2, 
                                                const model::MCInfo& info,
                                                bool is_excluded = false);

/**
 * @brief Apply periodic boundary conditions to distance vector
 * 
 * @param dx, dy, dz Distance components (modified in place)
 * @param box Box dimensions
 */
void applyPeriodicBoundaryConditions(float& dx, float& dy, float& dz, 
                                   const float box[3]);

/**
 * @brief Compute erfc-damped Coulomb interaction
 * 
 * @param r Distance between particles
 * @param q1, q2 Particle charges
 * @param alpha PME alpha parameter
 * @return Electrostatic energy contribution
 */
double computeErfcCoulomb(double r, double q1, double q2, double alpha);

/**
 * @brief Check if particle pair is within cutoff
 * 
 * @param r2 Squared distance
 * @param cutoff2 Squared cutoff distance
 * @return true if within cutoff
 */
inline bool withinCutoff(double r2, double cutoff2) {
    return r2 <= cutoff2;
}

/**
 * @brief Validate real space calculation parameters
 * 
 * @param params PME parameters
 * @return true if parameters are valid
 */
bool validateRealSpaceParameters(const PMEParams& params);

// <agent-hook:pme_realspace>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 