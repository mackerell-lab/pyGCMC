#ifndef PMEGRIDPREPARE_HPP
#define PMEGRIDPREPARE_HPP

#include "model/ModelModule.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Select atoms to process based on fixed_only flag
 * 
 * @param state MC state
 * @param fixed_only If true, only process fixed residues
 * @param atomsToProcess Output vector of atom indices
 * @param totalCharge Output total charge
 */
void selectAtomsToProcess(const model::MCState& state, bool fixed_only, 
                         std::vector<int>& atomsToProcess, double& totalCharge);

/**
 * @brief Calculate reciprocal lattice vectors
 * 
 * @param box Box dimensions
 * @param recipBoxVectors Output reciprocal box vectors [3][3]
 */
void calculateReciprocalLatticeVectors(const float* box, double recipBoxVectors[3][3]);

/**
 * @brief Output debug information about atoms and charges
 * 
 * @param state MC state
 * @param atomsToProcess List of atoms to process
 */
void outputDebugInfo(const model::MCState& state, const std::vector<int>& atomsToProcess);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEGRIDPREPARE_HPP 