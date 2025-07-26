#pragma once

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Complete PGP calculation including both electrostatic and LJ interactions
 * 
 * This function combines PGP electrostatic calculations with proper LJ handling,
 * similar to PME Complete but for the PGP method.
 */

/**
 * @brief Calculate complete PGP energy including intramolecular LJ
 * 
 * This function calculates:
 * 1. PGP electrostatic energy (grid interpolation + real space + self)
 * 2. Complete LJ interactions including intramolecular pairs
 * 
 * The key difference from regular PGP is that this includes ALL LJ interactions,
 * not just intermolecular ones, making it consistent with PME Complete.
 * 
 * @param state MC state containing system information
 */
void computeSystemEnergyPGPComplete(model::MCState& state);

/**
 * @brief Calculate movement energy using PGP Complete
 * 
 * This calculates energy for movement residues only, including:
 * 1. PGP electrostatic energy for movement residues
 * 2. LJ interactions between movement residues and all residues
 * 
 * @param state MC state containing system information
 */
void computeMovementEnergyPGPComplete(model::MCState& state);

/**
 * @brief Corrected implementation of movement energy using PGP Complete
 * 
 * This properly calculates energy for movement residues:
 * 1. Grid interpolation for movement atoms in fixed potential
 * 2. Real space interactions of movement atoms with ALL atoms
 * 3. LJ interactions of movement atoms with ALL atoms
 * 4. Self energy of movement atoms only
 * 
 * @param state MC state containing system information
 */
void computeMovementEnergyPGPCompleteCorrect(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc