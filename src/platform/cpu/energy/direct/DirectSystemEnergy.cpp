#include "DirectSystemEnergy.hpp"
#include "DirectCore.hpp"
#include "DirectPBCCalculation.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

void computeMovementEnergy(model::MCState& state) {
    computeNonbondedEnergy(state, false, true, false);  // No cutoff, movement residues only, no PBC
}

void computeMovementEnergyCutoff(model::MCState& state) {
    computeNonbondedEnergy(state, true, true, false);  // With cutoff, movement residues only, no PBC
}

void computeSystemEnergy(model::MCState& state) {
    computeNonbondedEnergy(state, false, false, false);  // No cutoff, all residues, no PBC
}

void computeSystemEnergyCutoff(model::MCState& state) {
    computeNonbondedEnergy(state, true, false, false);  // With cutoff, all residues, no PBC
}

void computeSystemVdwEnergyCutoff(model::MCState& state) {
    // Call computeNonbondedEnergy with VDW-only flag
    computeNonbondedEnergy(state, true, false, false, true);  // use_cutoff=true, movement_only=false, use_pbc=false, vdw_only=true
}

/**
 * @brief Unified system energy calculation function (Direct method)
 * 
 * Uses the direct calculation method to compute non-bonded interaction energies between all atoms in the system.
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    if (use_pbc) {
        if (use_cutoff) {
            computeSystemEnergyPBCCutoff(state);
        } else {
            computeSystemEnergyPBC(state);
        }
    } else {
        if (use_cutoff) {
            computeSystemEnergyCutoff(state);
        } else {
            // Use local function directly, no scope qualification needed
            // This avoids ambiguity issues
            computeNonbondedEnergy(state, false, false, false);
        }
    }
}

/**
 * @brief Unified energy calculation function for movement residues (Direct method)
 * 
 * Uses the direct calculation method to compute non-bonded interaction energies between movement residues and other atoms in the system.
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    if (use_pbc) {
        // Currently there is no dedicated PBC version of movement residue energy calculation function
        // We use the full system calculation, which might be slightly slower
        if (use_cutoff) {
            computeSystemEnergyPBCCutoff(state);
        } else {
            computeSystemEnergyPBC(state);
        }
    } else {
        if (use_cutoff) {
            computeMovementEnergyCutoff(state);
        } else {
            // Use local function directly, no scope qualification needed
            // This avoids ambiguity issues
            computeNonbondedEnergy(state, false, true, false);
        }
    }
}

/**
 * @brief Unified interface function for calculating van der Waals energy only
 * 
 * Uses the direct calculation method to compute only van der Waals interaction energies, without electrostatic energy.
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    if (use_cutoff) {
        computeSystemVdwEnergyCutoff(state);
    } else {
        // If there's no dedicated function for VDW energy calculation without cutoff, use the general function but only keep the VDW part
        computeNonbondedEnergy(state, use_cutoff, false, use_pbc, true);
    }
}

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 