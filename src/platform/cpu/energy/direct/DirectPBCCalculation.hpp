#pragma once

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

/**
 * @brief Calculate system energy using periodic boundary conditions (no cutoff)
 * 
 * @param state System state containing box dimensions and other parameters
 */
void computeSystemEnergyPBC(model::MCState& state);

/**
 * @brief Calculate system energy using periodic boundary conditions with cutoff
 * 
 * @param state System state containing box dimensions and cutoff parameters
 */
void computeSystemEnergyPBCCutoff(model::MCState& state);

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 