#include "PGPSystem.hpp"
#include "PGPGlobal.hpp"
#include "PGPCore.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "PGPInterpolation.hpp"
#include "PGPPrecompute.hpp"
#include "../lj/LJMain.hpp"
#include "../common/EnergyDirectCore.hpp"
#include "platform/Platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

static bool isMovementResidue(int residue_index, const model::MCState& state) {
    if (state.movementResidues.empty()) {
        return residue_index >= 0 &&
               residue_index < state.activeResidueCount &&
               !state.residues[residue_index].fixed;
    }

    for (const auto& movementInfo : state.movementResidues) {
        if (residue_index >= movementInfo.startIndex &&
            residue_index < movementInfo.startIndex + movementInfo.activeCount) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Calculate real-space PGP energy with correct erfc implementation
 * This fixes the bug where getPGPParams().erfcApprox uses pme_params tables
 */
static void computeRealSpacePGPFixed(model::MCState& state, bool movement_only, bool store_in_residues) {
    const auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float* box = state.info.box;
    const float cutoff = state.info.cutoff;
    const float cutoff2 = cutoff * cutoff;

    double real_space_total = 0.0;
    int debug_count = 0;
    const int max_debug_pairs = 5;

    // Calculate residue pairs.
    // For movement_only mode, evaluate pairs where at least one residue is movement.
    for(int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if (!residues[r1].active) continue;

        const bool r1_is_movement = movement_only ? isMovementResidue(r1, state) : true;
        if (movement_only && !r1_is_movement) continue;

        const int r2_start = movement_only ? 0 : r1;
        for(int r2 = r2_start; r2 < state.activeResidueCount; r2++) {
            if (!residues[r2].active) continue;

            const bool r2_is_movement = movement_only ? isMovementResidue(r2, state) : true;
            if (movement_only) {
                if (!r1_is_movement && !r2_is_movement) continue;
                // Moving-moving pairs are evaluated once.
                if (r2_is_movement && r2 < r1) continue;
            }

            // If both residues are fixed, skip (already in precomputed grid)
            if(residues[r1].fixed && residues[r2].fixed) continue;

            // Loop over atoms in each residue
            for(int i = residues[r1].atomStart;
                i < residues[r1].atomStart + residues[r1].atomCount; i++) {

                if(i >= state.activeAtomCount) continue;

                int j_start = (r1 == r2) ? i + 1 : residues[r2].atomStart;
                int j_end = residues[r2].atomStart + residues[r2].atomCount;

                for(int j = j_start; j < j_end; j++) {

                    if(j >= state.activeAtomCount) continue;

                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;

                    // Apply periodic boundary conditions
                    dx -= box[0] * round(dx / box[0]);
                    dy -= box[1] * round(dy / box[1]);
                    dz -= box[2] * round(dz / box[2]);

                    float r2 = dx*dx + dy*dy + dz*dz;

                    // Skip pairs beyond cutoff
                    if(r2 > cutoff2) continue;

                    // Compute energy
                    float r = sqrt(r2);
                    float qi = atoms[i].charge;
                    float qj = atoms[j].charge;

                    // Skip neutral atoms
                    if(std::abs(qi) < 1e-6 || std::abs(qj) < 1e-6) continue;

                    // FIXED: Calculate erfc directly using getPGPParams()
                    double alphaR = getPGPParams().alpha * r;
                    double erfc_val = std::erfc(alphaR);
                    double pair_energy = qi * qj * erfc_val / r;

                    // Print debug information
                    if (platform::is_debug_mode() && debug_count < max_debug_pairs) {
                        platform::log(LogLevel::DEBUG,
                            "PGPFixed: Atom pair (", i, ",", j, "): ",
                            "r = ", r, " nm, ",
                            "q1*q2 = ", qi * qj, ", ",
                            "erfc(", alphaR, ") = ", erfc_val, ", ",
                            "energy = ", pair_energy,
                            ", with COULOMB = ", COULOMB * pair_energy, " kJ/mol");
                        debug_count++;
                    }

                    // Accumulate to total energy
                    real_space_total += pair_energy;

                    // Store energy in residues
                    if (store_in_residues) {
                        if (r1 == r2) {
                            // Intra-residue: all energy goes to this residue
                            residues[r1].energy_elec += pair_energy;
                        } else {
                            // Inter-residue: split energy
                            residues[r1].energy_elec += pair_energy / 2.0f;
                            residues[r2].energy_elec += pair_energy / 2.0f;
                        }
                    }
                }
            }
        }
    }

    // Store total real-space energy (not yet multiplied by COULOMB)
    state.ewald_energy.real_space = real_space_total;
}

/**
 * @brief Fixed version of PGP system energy calculation
 * Uses correct erfc calculation and includes intra-residue interactions
 */
void computeSystemEnergyPGPFixed(model::MCState& state) {
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing total system energy using PGP-Fixed method");
    }

    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);

    // 2. Calculate real space part with FIXED implementation
    computeRealSpacePGPFixed(state, false, true);

    // 3. Calculate self energy correction
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, false);

    // 4. Calculate LJ interactions using direct cutoff method
    computeSystemVdwEnergyCutoff(state);

    // Multiply real space energy by COULOMB constant
    state.ewald_energy.real_space *= COULOMB;

    // Calculate total LJ energy
    double vdw_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            vdw_total += residue.energy_vdw;
        }
    }

    // Calculate total energy
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space +
                             state.ewald_energy.self + vdw_total;

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP-Fixed system energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

/**
 * @brief Fixed version of PGP movement energy calculation
 */
void computeMovementEnergyPGPFixed(model::MCState& state) {
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing movement residue energy using PGP-Fixed method");
    }

    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);

    // 2. Calculate real space part (only for moving residues)
    computeRealSpacePGPFixed(state, true, true);

    // 3. Calculate self energy correction (only for moving residues)
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, true);

    // 4. Calculate LJ interactions for movement residues only.
    // If movementResidues is empty, use all non-fixed active residues as implicit movement set.
    const bool injected_implicit_movement = state.movementResidues.empty();
    if (injected_implicit_movement) {
        state.movementResidues.clear();
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (!state.residues[i].active || state.residues[i].fixed) continue;
            model::MCMovementResidueInfo movement_info;
            movement_info.startIndex = i;
            movement_info.activeCount = 1;
            movement_info.totalCount = 1;
            movement_info.resName = "";
            state.movementResidues.push_back(movement_info);
        }
    }

    const int original_num_movement_types = state.forcefield.numMovementTypes;
    if (state.forcefield.numMovementTypes <= 0) {
        state.forcefield.numMovementTypes = state.forcefield.numTotalTypes;
    }
    computeMovementVdwEnergyDirect(state, true, false);
    state.forcefield.numMovementTypes = original_num_movement_types;

    // Multiply real space energy by COULOMB constant
    state.ewald_energy.real_space *= COULOMB;

    // Only accumulate LJ energy for moving residues
    double vdw_total = 0.0;
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                vdw_total += state.residues[i].energy_vdw;
            }
        }
    }

    // Calculate total energy
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space +
                             state.ewald_energy.self + vdw_total;

    if (injected_implicit_movement) {
        state.movementResidues.clear();
    }

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP-Fixed movement energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
