#include "DirectComposite.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

void DirectComposite::calculateSystemEnergy(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeSystemEnergyDirect(state, use_cutoff, use_pbc);
}

void DirectComposite::calculateMovementEnergy(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeMovementEnergyDirect(state, use_cutoff, use_pbc);
}

void DirectComposite::calculateVdwEnergy(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeSystemVdwEnergyDirect(state, use_cutoff, use_pbc);
}

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 