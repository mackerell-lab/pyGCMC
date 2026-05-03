#include "PGPSelf.hpp"
#include "PGPGlobal.hpp"
#include "PGPCore.hpp"
#include "platform/Platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate self energy correction for PGP method
 */
double computeSelfEnergyPGPImpl(model::MCState& state, bool movement_only) {
    double self_energy = 0.0;
    double sum_q2 = 0.0;

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing self energy for PGP method");
        platform::log(LogLevel::DEBUG, "Movement only: ", movement_only);
    }

    // Sum up squares of charges
    if (movement_only) {
        // Only include moving residues
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex;
                 i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if (!state.residues[i].active) continue;

                // Sum q^2 for all atoms in this movement residue
                for (int j = 0; j < state.residues[i].atomCount; j++) {
                    int atomIdx = state.residues[i].atomStart + j;
                    double q = state.atoms[atomIdx].charge;
                    sum_q2 += q * q;
                }
            }
        }
    }
    else {
        // Sum q^2 for active residues/atoms only.
        // Do not rely on activeAtomCount because deleted residues can leave ghost atoms.
        for (int r = 0; r < state.activeResidueCount; ++r) {
            const auto& residue = state.residues[r];
            if (!residue.active) continue;
            for (int j = residue.atomStart; j < residue.atomStart + residue.atomCount; ++j) {
                if (j < 0 || j >= static_cast<int>(state.atoms.size())) continue;
                double charge = state.atoms[j].charge;
                double q2 = charge * charge;
                sum_q2 += q2;
            }
        }
    }

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Sum of q² = ", sum_q2);
    }

    // Self-energy formula: -ONE_4PI_EPS0 * alpha / sqrt(PI) * sum_q2
    double prefactor = -COULOMB * getPGPParams().alpha / sqrt(M_PI);
    self_energy = prefactor * sum_q2;

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Self energy prefactor = ", prefactor,
                     ", resulting self energy = ", self_energy);
    }

    return self_energy;
}

/**
 * @brief Public interface wrapper for self energy PGP calculation
 */
double computeSelfEnergyPGP(model::MCState& state, bool movement_only) {
    return computeSelfEnergyPGPImpl(state, movement_only);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
