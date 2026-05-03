// PGPCompleteCorrected.cpp
// Correct implementation of PGP Complete based on PGP principles

#include "PGPComplete.hpp"
#include "PGPGlobal.hpp"
#include "PGPCore.hpp"
#include "PGPInterpolation.hpp"
#include "PGPSelf.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/Platform.hpp"
#include <cmath>
#include <cstdio>

namespace pygcmc {
namespace platform {
namespace cpu {

// Helper function to check if a residue is in movement list
static bool isMovementResidue(int res_idx, const model::MCState& state) {
    // If no movement residues specified, treat all non-fixed residues as movement
    if (state.movementResidues.empty()) {
        return !state.residues[res_idx].fixed;
    }

    for (const auto& movementInfo : state.movementResidues) {
        if (res_idx >= movementInfo.startIndex &&
            res_idx < movementInfo.startIndex + movementInfo.activeCount) {
            return true;
        }
    }
    return false;
}

// Calculate real space electrostatics following PGP principles
static void calculateRealSpacePGPComplete(model::MCState& state, bool movement_only) {
    const auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float* box = state.info.box;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const double alpha = getPGPParams().alpha;

    // Reset real space energy
    state.ewald_energy.real_space = 0.0;

    // Reset electrostatic energies in residues
    for (auto& res : residues) {
        if (res.active) {
            res.energy_elec = 0.0f;
        }
    }

    double real_space_total = 0.0;

    // Loop over residue pairs
    for (int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if (!residues[r1].active) continue;

        const bool r1_is_movement = movement_only ? isMovementResidue(r1, state) : true;
        if (movement_only && !r1_is_movement) continue;

        // In movement_only mode, movement-fixed interactions must include all r2.
        // Use r2_start=0 and keep unique counting for movement-movement pairs.
        const int r2_start = movement_only ? 0 : r1;
        for (int r2 = r2_start; r2 < state.activeResidueCount; r2++) {
            if (!residues[r2].active) continue;

            const bool r2_is_movement = movement_only ? isMovementResidue(r2, state) : true;
            if (movement_only) {
                if (!r1_is_movement && !r2_is_movement) continue;
                if (r2_is_movement && r2 < r1) continue;
            }

            // CRITICAL: If both residues are fixed, skip
            // Their interaction is already in the precomputed grid
            if (residues[r1].fixed && residues[r2].fixed) continue;

            // Loop over atom pairs
            for (int i = residues[r1].atomStart;
                 i < residues[r1].atomStart + residues[r1].atomCount; i++) {

                if (i >= state.activeAtomCount) continue;

                const double qi = atoms[i].charge;
                if (std::abs(qi) < 1e-6) continue;

                int j_start = (r1 == r2) ? i + 1 : residues[r2].atomStart;
                int j_end = residues[r2].atomStart + residues[r2].atomCount;

                for (int j = j_start; j < j_end; j++) {
                    if (j >= state.activeAtomCount) continue;

                    const double qj = atoms[j].charge;
                    if (std::abs(qj) < 1e-6) continue;

                    // Calculate distance with PBC
                    double dx = atoms[i].x - atoms[j].x;
                    double dy = atoms[i].y - atoms[j].y;
                    double dz = atoms[i].z - atoms[j].z;

                    // Apply minimum image convention
                    float dx_f = static_cast<float>(dx);
                    float dy_f = static_cast<float>(dy);
                    float dz_f = static_cast<float>(dz);
                    applyPBC(dx_f, dy_f, dz_f, box);
                    dx = dx_f; dy = dy_f; dz = dz_f;

                    const double r2_dist = dx*dx + dy*dy + dz*dz;

                    // Apply cutoff
                    if (r2_dist > cutoff2) continue;

                    // Skip extremely close atoms
                    if (r2_dist < 1e-12) continue;

                    const double r = std::sqrt(r2_dist);

                    // Calculate erfc(alpha*r)/r
                    const double alphar = alpha * r;
                    const double erfc_val = std::erfc(alphar);
                    const double pair_energy = qi * qj * erfc_val / r;

                    // Accumulate energy
                    real_space_total += pair_energy;

                    // Store in residues
                    if (r1 == r2) {
                        // Intra-residue: all energy to this residue
                        residues[r1].energy_elec += pair_energy;
                    } else {
                        // Inter-residue: split energy
                        residues[r1].energy_elec += pair_energy * 0.5;
                        residues[r2].energy_elec += pair_energy * 0.5;
                    }
                }
            }
        }
    }

    state.ewald_energy.real_space = real_space_total;
}

// Calculate LJ energy for PGP Complete
static void calculateLJPGPComplete(model::MCState& state, bool movement_only) {
    const auto& atoms = state.atoms;
    const auto& forcefield = state.forcefield;
    auto& residues = state.residues;
    const float* box = state.info.box;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const bool usePairtypes14 = forcefield.pairtypes14Enabled;

    platform::log(LogLevel::DEBUG, "calculateLJPGPComplete: movement_only=", movement_only);
    platform::log(LogLevel::DEBUG, "Active residue count: ", state.activeResidueCount);


    // Reset VDW energies
    for (auto& res : residues) {
        if (res.active) {
            res.energy_vdw = 0.0f;
        }
    }

    // Loop over residue pairs
    for (int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if (!residues[r1].active) continue;

        bool r1_is_movement = isMovementResidue(r1, state);

        // For movement_only mode, skip if r1 is not movement
        if (movement_only && !r1_is_movement) {
            platform::log(LogLevel::DEBUG, "Skipping residue ", r1, " - not movement");
            continue;
        }

        platform::log(LogLevel::DEBUG, "Processing residue r1=", r1);

        const int r2_start = movement_only ? 0 : r1;
        for (int r2 = r2_start; r2 < state.activeResidueCount; r2++) {
            if (!residues[r2].active) continue;

            const bool r2_is_movement = isMovementResidue(r2, state);
            if (movement_only) {
                if (!r1_is_movement && !r2_is_movement) continue;
                if (r2_is_movement && r2 < r1) continue;
            }

            platform::log(LogLevel::DEBUG, "Processing residue pair r1=", r1, " r2=", r2);

            // Loop over atom pairs
            for (int i = residues[r1].atomStart;
                 i < residues[r1].atomStart + residues[r1].atomCount; i++) {

                if (i >= state.activeAtomCount) continue;
                const int type_i = atoms[i].type;

                int j_start = (r1 == r2) ? i + 1 : residues[r2].atomStart;
                int j_end = residues[r2].atomStart + residues[r2].atomCount;

                for (int j = j_start; j < j_end; j++) {
                    if (j >= state.activeAtomCount) continue;
                    const int type_j = atoms[j].type;

                    // Get LJ parameters
                    int param_idx = type_i * forcefield.numTotalTypes + type_j;


                    if (param_idx >= static_cast<int>(forcefield.ljSigma.size())) {
                        continue;
                    }

                    double sigma = forcefield.ljSigma[param_idx];
                    double epsilon = forcefield.ljEps[param_idx];
                    if (usePairtypes14 && state.isPair14(i, j) &&
                        forcefield.hasPairtype14(type_i, type_j)) {
                        sigma = forcefield.ljSigma14[param_idx];
                        epsilon = forcefield.ljEps14[param_idx];
                    }


                    if (epsilon == 0.0 || sigma == 0.0) {
                        continue;
                    }

                    // Calculate distance with PBC
                    double dx = atoms[i].x - atoms[j].x;
                    double dy = atoms[i].y - atoms[j].y;
                    double dz = atoms[i].z - atoms[j].z;

                    // Apply minimum image convention
                    float dx_f = static_cast<float>(dx);
                    float dy_f = static_cast<float>(dy);
                    float dz_f = static_cast<float>(dz);
                    applyPBC(dx_f, dy_f, dz_f, box);
                    dx = dx_f; dy = dy_f; dz = dz_f;

                    const double r2_dist = dx*dx + dy*dy + dz*dz;

                    // Apply cutoff
                    if (r2_dist > cutoff2) continue;

                    // Skip extremely close atoms
                    if (r2_dist < 1e-12) continue;

                    // Calculate LJ energy
                    const double sigma2 = sigma * sigma;
                    const double sigma6 = sigma2 * sigma2 * sigma2;
                    const double sigma12 = sigma6 * sigma6;
                    const double r6 = r2_dist * r2_dist * r2_dist;
                    const double r12 = r6 * r6;

                    const double lj_energy = 4.0 * epsilon * (sigma12/r12 - sigma6/r6);

                    platform::log(LogLevel::DEBUG, "LJ pair i=", i, " j=", j, " r=", std::sqrt(r2_dist), " energy=", lj_energy);

                    // Store in residues
                    if (r1 == r2) {
                        // Intra-residue: all energy to this residue
                        residues[r1].energy_vdw += lj_energy;
                    } else {
                        // Inter-residue: split energy
                        residues[r1].energy_vdw += lj_energy * 0.5;
                        residues[r2].energy_vdw += lj_energy * 0.5;
                    }
                }
            }
        }
    }
}

// Corrected implementation of movement energy using PGP Complete
void computeMovementEnergyPGPCompleteCorrect(model::MCState& state) {
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }


    platform::log(LogLevel::DEBUG, "Computing movement energy using corrected PGP Complete method");
    platform::log(LogLevel::DEBUG, "Number of movement residue blocks: ", state.movementResidues.size());

    // If no movement residues specified, treat all non-fixed residues as movement
    if (state.movementResidues.empty()) {
        platform::log(LogLevel::DEBUG, "No movement residues specified, treating all non-fixed as movement");
    }

    // 1. Calculate grid potential interpolation
    // This gives the interaction of movement atoms with the fixed potential
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);

    // 2. Calculate real space electrostatics
    // For movement residues interacting with all residues
    // Skip fixed-fixed pairs (already in grid)
    calculateRealSpacePGPComplete(state, true);

    // 3. Calculate self energy for movement atoms only
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, true);

    // 4. Calculate LJ for movement residues
    calculateLJPGPComplete(state, true);

    // Apply COULOMB constant to real space and residue energies
    state.ewald_energy.real_space *= COULOMB;
    for (auto& res : state.residues) {
        if (res.active) {
            res.energy_elec *= COULOMB;
        }
    }

    // Calculate total VDW from movement residues only
    double vdw_total = 0.0;

    if (state.movementResidues.empty()) {
        // If no movement residues specified, sum VDW from all non-fixed residues
        for (int i = 0; i < state.activeResidueCount; i++) {
            if (state.residues[i].active && !state.residues[i].fixed) {
                vdw_total += state.residues[i].energy_vdw;
            }
        }
    } else {
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex;
                 i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if (state.residues[i].active) {
                    vdw_total += state.residues[i].energy_vdw;
                }
            }
        }
    }


    // Set reciprocal (grid) and calculate total
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space +
                             state.ewald_energy.self + vdw_total;

    platform::log(LogLevel::INFO, "PGP Complete Movement (corrected v2):");
    platform::log(LogLevel::INFO, "  Grid (reciprocal): ", grid_energy);
    platform::log(LogLevel::INFO, "  Real space: ", state.ewald_energy.real_space);
    platform::log(LogLevel::INFO, "  Self: ", state.ewald_energy.self);
    platform::log(LogLevel::INFO, "  VDW: ", vdw_total);
    platform::log(LogLevel::INFO, "  Total: ", state.ewald_energy.total);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
