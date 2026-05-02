#include "PGPReal.hpp"
#include "PGPGlobal.hpp"
#include "PGPCore.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate real-space part of the PGP method for short-range electrostatics
 */
void computeRealSpacePGPImpl(model::MCState& state, bool movement_only, bool store_in_residues) {
    const auto& box = state.info.box;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float cutoff2 = getPGPParams().cutoff * getPGPParams().cutoff;

    // Reset electrostatic energy
    for(auto& residue : residues) {
        if(residue.active) {
            residue.energy_elec = 0.0f;
        }
    }

    // Real-space total energy
    double real_space_total = 0.0;

    // Add debug information
    int debug_count = 0;
    const int max_debug_pairs = 5;

    // Loop over all residue pairs
    for(int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if(!residues[r1].active) continue;
        if(movement_only) {
            bool in_movement = false;
            for(const auto& movementInfo : state.movementResidues) {
                if(r1 >= movementInfo.startIndex &&
                   r1 < movementInfo.startIndex + movementInfo.activeCount) {
                    in_movement = true;
                    break;
                }
            }
            if(!in_movement) continue;
        }

        for(int r2 = r1 + 1; r2 < state.activeResidueCount; r2++) {
            if(!residues[r2].active) continue;

            // If both residues are fixed, their interaction is already included in the precomputed grid potential
            if(residues[r1].fixed && residues[r2].fixed) continue;

            // Loop over atoms in each residue
            for(int i = residues[r1].atomStart;
                i < residues[r1].atomStart + residues[r1].atomCount; i++) {
                // Ensure atom index is valid
                if(i >= state.activeAtomCount) continue;

                for(int j = residues[r2].atomStart;
                    j < residues[r2].atomStart + residues[r2].atomCount; j++) {
                    // Ensure atom index is valid
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

                    // Calculate real space contribution for PGP - only erfc part
                    double term = getPGPParams().erfcApprox(r);
                    double pair_energy = qi * qj * term / r;

                    // Print debug information
                    if (platform::is_debug_mode() && debug_count < max_debug_pairs) {
                        platform::log(LogLevel::DEBUG,
                            "Debug energyPGP: Atom pair (", i, ",", j, "): ",
                            "r = ", r, " nm, ",
                            "q1*q2 = ", qi * qj, ", ",
                            "erfc term = ", term, ", ",
                            "energy = ", pair_energy,
                            ", with COULOMB = ", COULOMB * pair_energy, " kJ/mol");
                        debug_count++;
                    }

                    // Accumulate to total energy
                    real_space_total += pair_energy;

                    // Decide how to store energy based on parameters
                    if (store_in_residues) {
                        // Each residue gets half of the pair interaction energy
                        residues[r1].energy_elec += pair_energy / 2.0f;
                        residues[r2].energy_elec += pair_energy / 2.0f;
                    }
                }
            }
        }
    }

    // Store total real-space energy (not yet multiplied by COULOMB)
    state.ewald_energy.real_space = real_space_total;
}

/**
 * @brief Public interface wrapper for real-space PGP calculation
 */
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues) {
    computeRealSpacePGPImpl(state, movement_only, store_in_residues);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
