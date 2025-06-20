#ifndef PGPSELFENERGY_HPP
#define PGPSELFENERGY_HPP

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate self energy correction for PGP method
 * 
 * @param state MC state
 * @param movement_only Only calculate for moving residues
 * @return Self energy correction
 */
double computeSelfEnergyPGPImpl(model::MCState& state, bool movement_only);

/**
 * @brief Public interface wrapper for self energy PGP calculation
 */
double computeSelfEnergyPGP(model::MCState& state, bool movement_only);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PGPSELFENERGY_HPP 