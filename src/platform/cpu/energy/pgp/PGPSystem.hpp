#ifndef PGPSYSTEMENERGY_HPP
#define PGPSYSTEMENERGY_HPP

#include "model/ModelModule.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Function Location Guide for AI Agents:
 * - System energy calculation: computeSystemEnergyPGP*
 * - Movement energy calculation: computeMovementEnergyPGP*
 * - Grid interpolation: PGPInterpolation.hpp -> interpolateMoleculeEnergy*, calculateMoleculeEnergy*
 * - Global energy calculation: PGPInterpolation.hpp -> computeMoleculeEnergyGlobal*
 */

// === System Energy Calculation Functions ===

/**
 * @brief Use PGP method to calculate system energy
 * 
 * @param state MC state
 */
void computeSystemEnergyPGPImpl(model::MCState& state);

/**
 * @brief Use PGP method to calculate energy of moving residues
 * 
 * @param state MC state
 */
void computeMovementEnergyPGPImpl(model::MCState& state);

/**
 * @brief Public interface wrapper for system energy calculation
 */
void computeSystemEnergyPGP(model::MCState& state);

/**
 * @brief Public interface wrapper for movement energy calculation
 */
void computeMovementEnergyPGP(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PGPSYSTEMENERGY_HPP 