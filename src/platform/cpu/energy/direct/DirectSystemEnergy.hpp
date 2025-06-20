#pragma once

#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

// Basic direct calculation functions (following existing naming pattern)
void computeMovementEnergy(model::MCState& state);
void computeMovementEnergyCutoff(model::MCState& state);
void computeSystemEnergy(model::MCState& state);
void computeSystemEnergyCutoff(model::MCState& state);
void computeSystemVdwEnergyCutoff(model::MCState& state);

// Unified Direct method interfaces
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc);

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 