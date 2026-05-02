#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_MOVE_COMMON_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_MOVE_COMMON_HPP

#include "model/montecarlo/MCMain.hpp"
#include "MovementUtils.hpp"
#include <algorithm>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace move_common {

inline int countActiveResiduesOfType(
    const model::montecarlo::MCState& state,
    int moleculeType) {
    if (moleculeType < 0) {
        return state.activeResidueCount;
    }

    int count = 0;
    const int maxResidues = std::min(
        state.activeResidueCount,
        static_cast<int>(state.residues.size()));
    for (int i = 0; i < maxResidues; ++i) {
        const auto& residue = state.residues[i];
        if (residue.active && residue.type == moleculeType) {
            ++count;
        }
    }
    return count;
}

inline double sumResiduePairEnergy(
    const model::montecarlo::MCState& state,
    int excludedResidue = -1,
    bool activeOnly = false) {
    double total = 0.0;
    const int maxResidues = std::min(
        state.activeResidueCount,
        static_cast<int>(state.residues.size()));
    for (int i = 0; i < maxResidues; ++i) {
        if (i == excludedResidue) {
            continue;
        }
        const auto& residue = state.residues[i];
        if (activeOnly && !residue.active) {
            continue;
        }
        total += residue.energy_vdw;
        total += residue.energy_elec;
    }
    return 0.5 * total;
}

inline int selectRandomActiveResidue(const model::montecarlo::MCState& state) {
    const int maxResidues = std::min(
        state.activeResidueCount,
        static_cast<int>(state.residues.size()));
    if (maxResidues <= 0) {
        return -1;
    }

    std::vector<int> activeResidues;
    activeResidues.reserve(maxResidues);
    for (int i = 0; i < maxResidues; ++i) {
        if (state.residues[i].active) {
            activeResidues.push_back(i);
        }
    }

    if (activeResidues.empty()) {
        return -1;
    }

    const int idx = utils::RandomUtils::uniformInt(
        0,
        static_cast<int>(activeResidues.size()) - 1);
    return activeResidues[idx];
}

} // namespace move_common
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_MOVE_COMMON_HPP
