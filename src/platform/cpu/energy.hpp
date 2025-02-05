// src/platform/cpu/energy.hpp

#pragma once

#include "../../model/montecarlo.hpp"

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
 * @return The total energy (sum of all residue energies divided by 2).
 */
float computeNaiveNonbondedEnergy(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc