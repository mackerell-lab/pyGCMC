#pragma once

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Function Location Guide for AI Agents:
 * - System energy calculation: computeSystemEnergyPME
 * - Movement energy calculation: computeMovementEnergyPME  
 * - Component calculations: computeReciprocalPME, computeSelfEnergyPME, computeRealSpacePME
 * - Pair energy calculation: calcPairEnergyPME
 */

// === High-Level Energy Calculation Interfaces ===

/**
 * @brief Compute total system energy using PME
 * 
 * Calculates all components: real space + reciprocal space + self energy + VdW
 * 
 * @param state MC state containing system information
 */
void computeSystemEnergyPME(model::MCState& state);

/**
 * @brief Compute energy for moving residues only
 * 
 * Optimized calculation for Monte Carlo moves that only affect
 * a subset of the system.
 * 
 * @param state MC state
 */
void computeMovementEnergyPME(model::MCState& state);

// === Component Energy Calculation Functions ===

/**
 * @brief Compute reciprocal space energy using PME
 * 
 * @param state MC state
 * @return Reciprocal space energy
 */
double computeReciprocalPME(model::MCState& state);

/**
 * @brief Compute self-energy correction for PME
 * 
 * @param state MC state containing particle charges
 * @param movement_only If true, only calculate for moving particles
 * @return Self-energy correction (negative value)
 */
double computeSelfEnergyPME(model::MCState& state, bool movement_only);

/**
 * @brief Calculate real-space part of PME
 * 
 * Computes the real-space energy contribution using erfc(αr)/r terms
 * for all particle pairs within the cutoff distance.
 * 
 * @param state MC state containing system information
 * @param movement_only Whether to calculate only for moving residues
 * @param store_in_residues Whether to store energy in residue objects
 */
void computeRealSpacePME(model::MCState& state, bool movement_only, bool store_in_residues = true);

// === Pair Energy Calculation ===

/**
 * @brief Calculate pair energy for PME real-space part
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
std::pair<double, double> calcPairEnergyPME(double r2, double sigma, double eps, 
                                          double q1, double q2, 
                                          const model::MCInfo& info, 
                                          bool is_excluded = false);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 