#ifndef PGPSYSTEMENERGY_HPP
#define PGPSYSTEMENERGY_HPP

#include "model/montecarlo.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Function Location Guide for AI Agents:
 * - System energy calculation: computeSystemEnergyPGP*
 * - Movement energy calculation: computeMovementEnergyPGP*
 * - Grid interpolation: interpolateMoleculeEnergy*, calculateMoleculeEnergy*
 * - Global energy calculation: computeMoleculeEnergyGlobal*
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

// === Grid Interpolation Functions ===

/**
 * @brief Calculate moving molecule energy through interpolation
 * 
 * This function uses B-spline interpolation from the precomputed grid potential
 * to quickly evaluate the energy of moving molecules, avoiding direct calculation
 * of intermolecular interactions.
 * 
 * @param state MC state containing moving molecules and precomputed potential grid
 * @param energy Output parameter that stores the calculated energy value
 */
void interpolateMoleculeEnergyImpl(model::MCState& state, double& energy);

/**
 * @brief Calculate moving molecule energy through interpolation and return result
 * 
 * Wrapper function for interpolateMoleculeEnergyImpl that directly returns 
 * the calculated energy value. Convenient for Python calls and testing.
 * 
 * @param state MC state containing moving molecules and precomputed potential grid
 * @return Calculated energy value
 */
double calculateMoleculeEnergyImpl(model::MCState& state);

/**
 * @brief Global energy calculation with specified residue lists
 * 
 * This function calculates molecule energy for specified movement residues
 * and considers nearby residues, with multi-threading support.
 * 
 * @param state MC state
 * @param movementResidues List of residue indices to calculate energy for
 * @param nearbyResidues List of nearby residue indices for corrections
 * @param threadIndex Thread index for multi-threaded calculations
 * @return Calculated energy value
 */
double computeMoleculeEnergyGlobalImpl(model::MCState& state, 
                                       const std::vector<int>& movementResidues, 
                                       const std::vector<int>& nearbyResidues, 
                                       int threadIndex);

// === Public Interface Wrappers ===

/**
 * @brief Public interface for molecule energy interpolation
 */
void interpolateMoleculeEnergy(model::MCState& state, double& energy);

/**
 * @brief Public interface for molecule energy calculation
 */
double calculateMoleculeEnergy(model::MCState& state);

/**
 * @brief Public interface for global molecule energy calculation
 */
double computeMoleculeEnergyGlobal(model::MCState& state, 
                                   const std::vector<int>& movementResidues, 
                                   const std::vector<int>& nearbyResidues, 
                                   int threadIndex);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PGPSYSTEMENERGY_HPP 