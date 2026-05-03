#include "PMEGridPrep.hpp"
#include "PMEGlobal.hpp"
#include "PMECore.hpp"
#include "platform/Platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Select atoms to process based on fixed_only flag
 */
void selectAtomsToProcess(const model::MCState& state, bool fixed_only,
                         std::vector<int>& atomsToProcess, double& totalCharge) {
    const auto& atoms = state.atoms;
    const auto& residues = state.residues;

    atomsToProcess.clear();
    atomsToProcess.reserve(state.activeAtomCount);
    totalCharge = 0.0;

    if (fixed_only) {
        // Only include atoms from fixed residues
        for (int i = 0; i < state.activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (residue.active && residue.fixed) {
                // Add all atoms in this fixed residue
                for (int j = 0; j < residue.atomCount; ++j) {
                    int atomIndex = residue.atomStart + j;
                    atomsToProcess.push_back(atomIndex);
                    totalCharge += atoms[atomIndex].charge;
                }
            }
        }

        platform::log(LogLevel::INFO, "Processing ", atomsToProcess.size(),
                     " atoms from fixed residues, total charge: ", totalCharge);
    } else {
        // Process atoms from all active residues.
        // NOTE: gcmc_cpu keeps deleted atoms in the global arrays and does not always decrement
        // activeAtomCount; iterating by residues avoids accidentally including "ghost" atoms.
        for (int i = 0; i < state.activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (!residue.active) continue;
            for (int j = 0; j < residue.atomCount; ++j) {
                const int atomIndex = residue.atomStart + j;
                if (atomIndex < 0 || atomIndex >= static_cast<int>(atoms.size())) {
                    continue;
                }
                atomsToProcess.push_back(atomIndex);
                totalCharge += atoms[atomIndex].charge;
            }
        }

        platform::log(LogLevel::INFO, "Total system charge: ", totalCharge);
        platform::log(LogLevel::DEBUG, "Total system charge: " + std::to_string(totalCharge));
    }
}

/**
 * @brief Calculate reciprocal lattice vectors
 */
void calculateReciprocalLatticeVectors(const float* box, double recipBoxVectors[3][3]) {
    // Handle periodic box vectors - simplify for diagonal boxes
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };

    // Check if it's a diagonal box - no need to output
    bool isDiagonalBox = true;  // Always true since we force it to be diagonal

    // Initialize reciprocal vectors
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            recipBoxVectors[i][j] = 0.0;
        }
    }

    if (isDiagonalBox) {
        // Simple calculations for diagonal boxes - consistent with OpenMM (no 2π factor)
        recipBoxVectors[0][0] = 1.0 / box[0]; // 1/a
        recipBoxVectors[1][1] = 1.0 / box[1]; // 1/b
        recipBoxVectors[2][2] = 1.0 / box[2]; // 1/c
    } else {
        // Non-diagonal boxes require full calculation of reciprocal lattice vectors
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);

        // Calculate cross products and reciprocal lattice vectors - consistent with OpenMM (no 2π factor)
        recipBoxVectors[0][0] = (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[0][1] = (periodicBoxVectors[0][2] * periodicBoxVectors[2][1] - periodicBoxVectors[0][1] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[0][2] = (periodicBoxVectors[0][1] * periodicBoxVectors[1][2] - periodicBoxVectors[0][2] * periodicBoxVectors[1][1]) / det;

        recipBoxVectors[1][0] = (periodicBoxVectors[1][2] * periodicBoxVectors[2][0] - periodicBoxVectors[1][0] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[1][1] = (periodicBoxVectors[0][0] * periodicBoxVectors[2][2] - periodicBoxVectors[0][2] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[1][2] = (periodicBoxVectors[0][2] * periodicBoxVectors[1][0] - periodicBoxVectors[0][0] * periodicBoxVectors[1][2]) / det;

        recipBoxVectors[2][0] = (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[2][1] = (periodicBoxVectors[0][1] * periodicBoxVectors[2][0] - periodicBoxVectors[0][0] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[2][2] = (periodicBoxVectors[0][0] * periodicBoxVectors[1][1] - periodicBoxVectors[0][1] * periodicBoxVectors[1][0]) / det;
    }
}

/**
 * @brief Output debug information about atoms and charges
 */
void outputDebugInfo(const model::MCState& state, const std::vector<int>& atomsToProcess) {
    if (!platform::is_debug_mode()) {
        return;
    }

    platform::log(LogLevel::DEBUG, "Spreading charges onto PME grid");

    // Calculate total system charge
    double totalCharge = 0.0;
    for (int i = 0; i < state.activeAtomCount; i++) {
        totalCharge += state.atoms[i].charge;
    }
    platform::log(LogLevel::DEBUG, "Total system charge: " + std::to_string(totalCharge));

    // Output initial grid values
    platform::log(LogLevel::DEBUG, "Initial values of the first 10 grid points:");
    for (int i = 0; i < 10 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
        platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(i) + "] = " + std::to_string(pme_params.pmeGrid[i].real()));
    }

    // Output charge values
    platform::log(LogLevel::DEBUG, "Charge values of the first 10 atoms:");
    for (int i = 0; i < 10 && i < state.activeAtomCount; i++) {
        platform::log(LogLevel::DEBUG, "  Atom[" + std::to_string(i) + "] charge = " + std::to_string(state.atoms[i].charge));
    }

    // Print some atom charge values to verify if there are non-zero charges
    int numToPrint = std::min(10, static_cast<int>(atomsToProcess.size()));
    platform::log(LogLevel::DEBUG, "Charge values of the first " + std::to_string(numToPrint) + " processed atoms:");
    for (int i = 0; i < numToPrint; i++) {
        int atomIndex = atomsToProcess[i];
        platform::log(LogLevel::DEBUG, "  Atom[" + std::to_string(atomIndex) + "] charge = " + std::to_string(state.atoms[atomIndex].charge));
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
