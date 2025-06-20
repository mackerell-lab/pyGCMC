#include "DirectCore.hpp"
#include "DirectResidueEnergy.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <stdexcept>
#include <sstream>
#include <limits>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

/**
 * @brief Universal function for calculating all nonbonded interactions
 */
void computeNonbondedEnergy(model::MCState& state, bool use_cutoff, bool movement_only, bool use_pbc, bool vdw_only) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting nonbonded energy calculation ===";
        ss << "\nSystem state info:";
        ss << "\n  Active residue count: " << state.activeResidueCount;
        if (use_cutoff) {
            ss << "\n  Cutoff distance: " << state.info.cutoff << " nm";
        }
        if (movement_only) {
            ss << "\n  Calculating only for movement residues";
        }
        if (use_pbc) {
            ss << "\n  Using periodic boundary conditions";
        }
        if (vdw_only) {
            ss << "\n  Calculating VDW interactions only";
        }
        platform::log(LogLevel::DEBUG, ss.str());
    }

    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;

    // Validate basic state parameters
    if (state.activeResidueCount < 0 || 
        static_cast<size_t>(state.activeResidueCount) > residues.size()) {
        throw std::runtime_error("Invalid activeResidueCount: " + 
                               std::to_string(state.activeResidueCount) +
                               " (residues size: " + std::to_string(residues.size()) + ")");
    }

    if (forcefield.numTotalTypes <= 0) {
        throw std::runtime_error("Invalid numTotalTypes: " + 
                               std::to_string(forcefield.numTotalTypes));
    }

    if (movement_only && forcefield.numMovementTypes <= 0) {
        throw std::runtime_error("Invalid numMovementTypes: " + 
                               std::to_string(forcefield.numMovementTypes));
    }

    // Always expect full matrix size
    size_t expected_size = static_cast<size_t>(forcefield.numTotalTypes) * 
                           static_cast<size_t>(forcefield.numTotalTypes);
    
    if (forcefield.ljEps.size() != expected_size || forcefield.ljSigma.size() != expected_size) {
        std::stringstream ss;
        ss << "Force field parameters array size mismatch. Expected size "
           << expected_size
           << " (numTotalTypes * numTotalTypes), but got eps=" << forcefield.ljEps.size()
           << " sigma=" << forcefield.ljSigma.size();
        throw std::runtime_error(ss.str());
    }

    // Reset energies for all residues
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        if (!vdw_only) {  // Only reset electrostatic energy if we're calculating it
            residue.energy_elec = 0.0f;
        }
    }

    if (movement_only) {
        // Validate movement residues
        if (state.movementResidues.empty()) {
            throw std::runtime_error("No movement residues defined");
        }

        // Calculate energies only for movement residues
        for (const auto& movementInfo : state.movementResidues) {
            if (getEnergyDebugOutput()) {
                platform::log(LogLevel::DEBUG, "\nProcessing movement residue group: ", movementInfo.resName);
                platform::log(LogLevel::DEBUG, "  Start index: ", movementInfo.startIndex);
                platform::log(LogLevel::DEBUG, "  Active count: ", movementInfo.activeCount);
                platform::log(LogLevel::DEBUG, "  Total count: ", movementInfo.totalCount);
            }

            // Validate movement residue indices
            if (movementInfo.startIndex < 0 || 
                movementInfo.startIndex + movementInfo.activeCount > state.activeResidueCount) {
                throw std::runtime_error("Invalid movement residue range: [" + 
                                       std::to_string(movementInfo.startIndex) + ", " +
                                       std::to_string(movementInfo.startIndex + movementInfo.activeCount) + 
                                       ") exceeds active residue count " +
                                       std::to_string(state.activeResidueCount));
            }

            // Process active movement residues
            for (int i = movementInfo.startIndex;
                 i < movementInfo.startIndex + movementInfo.activeCount;
                 ++i) {
                if (!residues[i].active) continue;

                if (getEnergyDebugOutput()) {
                    platform::log(LogLevel::DEBUG, "\nProcessing movement residue ", i);
                    platform::log(LogLevel::DEBUG, "  Atom start: ", residues[i].atomStart);
                    platform::log(LogLevel::DEBUG, "  Atom count: ", residues[i].atomCount);
                }

                // Validate atom indices
                if (residues[i].atomStart < 0 || 
                    static_cast<size_t>(residues[i].atomStart + residues[i].atomCount) > atoms.size()) {
                    throw std::runtime_error("Invalid atom range for residue " + 
                                           std::to_string(i) + ": [" +
                                           std::to_string(residues[i].atomStart) + ", " +
                                           std::to_string(residues[i].atomStart + residues[i].atomCount) + 
                                           ") exceeds atoms size " +
                                           std::to_string(atoms.size()));
                }

                computeResidueNonbondedEnergy(state, i, use_cutoff, use_pbc);
            }
        }
    } else {
        // Calculate energies for all active residues
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (!residues[i].active) continue;

            // Validate atom indices
            if (residues[i].atomStart < 0 || 
                static_cast<size_t>(residues[i].atomStart + residues[i].atomCount) > atoms.size()) {
                throw std::runtime_error("Invalid atom range for residue " + 
                                       std::to_string(i) + ": [" +
                                       std::to_string(residues[i].atomStart) + ", " +
                                       std::to_string(residues[i].atomStart + residues[i].atomCount) + 
                                       ") exceeds atoms size " +
                                       std::to_string(atoms.size()));
            }

            computeResidueNonbondedEnergy(state, i, use_cutoff, use_pbc);
        }
    }

    // Output final energies if debug is enabled
    if (getEnergyDebugOutput()) {
        platform::log(LogLevel::DEBUG, "\n=== Final energies for all residues ===");
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (residues[i].active) {
                total_vdw += residues[i].energy_vdw;
                total_elec += residues[i].energy_elec;
                platform::log(LogLevel::DEBUG, "Residue ", i, ":");
                platform::log(LogLevel::DEBUG, "  VDW energy: ", residues[i].energy_vdw, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Electrostatic energy: ", residues[i].energy_elec, " kJ/mol");
                platform::log(LogLevel::DEBUG, "  Total energy: ", 
                            (residues[i].energy_vdw + residues[i].energy_elec), " kJ/mol");
            }
        }
        platform::log(LogLevel::DEBUG, "\nTotal system energy:");
        platform::log(LogLevel::DEBUG, "  VDW: ", total_vdw, " kJ/mol");
        platform::log(LogLevel::DEBUG, "  Electrostatic: ", total_elec, " kJ/mol");
        platform::log(LogLevel::DEBUG, "  Total: ", (total_vdw + total_elec), " kJ/mol");
        platform::log(LogLevel::DEBUG, "\n=== Completed nonbonded energy calculation ===");
    }
}

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 