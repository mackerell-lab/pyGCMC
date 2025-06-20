#ifndef PGPREALSPACE_HPP
#define PGPREALSPACE_HPP

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate real-space part of the PGP method for short-range electrostatics
 * 
 * @param state MC state
 * @param movement_only Only calculate for moving residues
 * @param store_in_residues Whether to store energy in residue objects
 */
void computeRealSpacePGPImpl(model::MCState& state, bool movement_only, bool store_in_residues);

/**
 * @brief Public interface wrapper for real-space PGP calculation
 */
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PGPREALSPACE_HPP 