// src/platform/cpu/energy.hpp

#pragma once

#include "../../model/montecarlo.hpp"
#include "../platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute a naive (simplified) nonbonded energy based on Lennard-Jones potential,
 *        considering only active movement molecules and all other active molecules.
 *
 * This implementation uses a simplified approach: 
 * - It computes actual distances between atoms
 * - For each atom pair, computes both vdw and electrostatic energies:
 *   - vdw: V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
 *   - elec: V = q1*q2/r
 * - Accumulates energies into each residue's energy_vdw and energy_elec parameters
 * 
 * @param state The current MC state containing residues and force field parameters.
 *              The residues' energy parameters will be modified.
 */
void computeNaiveNonbondedEnergy(model::MCState& state);

/**
 * @brief Compute nonbonded energy for all active residues in the system
 * 
 * This implementation:
 * - Computes pairwise interactions between all active residues
 * - For each atom pair, computes both vdw and electrostatic energies:
 *   - vdw: V = 4ε[(σ/r)¹² - (σ/r)⁶]
 *   - elec: V = k_c * q1*q2/r
 * - Accumulates energies into each residue's energy_vdw and energy_elec parameters
 * - Uses full force field parameter matrix (numTotalTypes × numTotalTypes)
 * - Implements safety checks for minimum distance and maximum energy
 * 
 * @param state The current MC state containing residues and force field parameters.
 *              The residues' energy parameters will be modified.
 */
void computeAllNonbondedEnergy(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc