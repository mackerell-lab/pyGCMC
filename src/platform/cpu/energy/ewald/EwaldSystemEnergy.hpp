#ifndef EWALDSYSTEMENERGY_HPP
#define EWALDSYSTEMENERGY_HPP

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute total system energy using Ewald summation
 * 
 * @param state MC state
 */
void computeSystemEnergy(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // EWALDSYSTEMENERGY_HPP 