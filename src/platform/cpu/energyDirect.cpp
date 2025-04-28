// src/platform/cpu/energyDirect.cpp

#include "energyDirect.hpp"
#include "energyLJ.hpp"  // 添加对新文件的引用
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>  // For output formatting

namespace pygcmc {
namespace platform {
namespace cpu {

// Now using the common debug flag instead of a local variable
// static bool debug_output = false;

/**
 * @brief Calculate LJ and Coulomb energy with safety checks
 */
std::pair<double, double> calcPairEnergy(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool calc_coulomb) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\n=== calcPairEnergy called ===";
        ss << "\nInput parameters:"
           << "\n  Distance² = " << r2 << " nm²"
           << "\n  Sigma = " << sigma << " nm"
           << "\n  Epsilon = " << eps << " kJ/mol"
           << "\n  q1 = " << q1 << " e"
           << "\n  q2 = " << q2 << " e";
        platform::log(LogLevel::DEBUG, ss.str());
    }

    // 使用简化版本的调用
    double vdw_energy = calcLJEnergy(r2, sigma, eps, info);
    
    double r = std::sqrt(r2);
    
    if (getEnergyDebugOutput()) {
        platform::log(LogLevel::DEBUG, "Distance r = ", r, " nm");
    }
    
    // Calculate Coulomb energy only if requested
    double elec_energy = 0.0;
    if (calc_coulomb) {
        elec_energy = COULOMB * q1 * q2 / r;  // kJ/mol
    }

    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\nEnergy calculation details:";
        ss << "\n  COULOMB constant = " << COULOMB;
        ss << "\n  q1*q2 = " << (q1 * q2);
        ss << "\nInitial energies:";
        ss << "\n  VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Electrostatic energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Apply energy capping for numerical stability
    if (getEnergyDebugOutput() && (std::abs(vdw_energy) > MAX_SAFE_ENERGY || std::abs(elec_energy) > MAX_SAFE_ENERGY)) {
        std::stringstream ss;
        ss << "\nEnergy capping applied:";
        ss << "\n  Original VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Original Elec energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // LJ能量已经在calculateLJEnergy中处理过限制，这里只需要处理elec_energy
    elec_energy = std::min(elec_energy, static_cast<double>(MAX_SAFE_ENERGY));
    elec_energy = std::max(elec_energy, -static_cast<double>(MAX_SAFE_ENERGY));
    
    if (getEnergyDebugOutput() && (std::abs(vdw_energy) > MAX_SAFE_ENERGY || std::abs(elec_energy) > MAX_SAFE_ENERGY)) {
        std::stringstream ss;
        ss << "\nAfter individual capping:";
        ss << "\n  Capped VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Capped Elec energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Cap total energy
    double total_energy = vdw_energy + elec_energy;
    double original_total = total_energy;
    float max_safe = MAX_SAFE_ENERGY;  // 使用float而不转换为double
    
    if (total_energy > max_safe) {
        // 使用浮点版本的safe值，减少转换
        double scale = max_safe / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (getEnergyDebugOutput()) {
            std::stringstream ss;
            ss << "\nTotal energy exceeded MAX_SAFE_ENERGY:";
            ss << "\n  Original total = " << original_total << " kJ/mol";
            ss << "\n  Scale factor = " << scale;
            ss << "\n  Final VDW = " << vdw_energy << " kJ/mol";
            ss << "\n  Final Elec = " << elec_energy << " kJ/mol";
            ss << "\n  Final total = " << (vdw_energy + elec_energy) << " kJ/mol";
            platform::log(LogLevel::DEBUG, ss.str());
        }
    } else if (total_energy < -max_safe) {
        // 使用浮点版本的safe值，减少转换
        double scale = -max_safe / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (getEnergyDebugOutput()) {
            std::stringstream ss;
            ss << "\nTotal energy below -MAX_SAFE_ENERGY:";
            ss << "\n  Original total = " << original_total << " kJ/mol";
            ss << "\n  Scale factor = " << scale;
            ss << "\n  Final VDW = " << vdw_energy << " kJ/mol";
            ss << "\n  Final Elec = " << elec_energy << " kJ/mol";
            ss << "\n  Final total = " << (vdw_energy + elec_energy) << " kJ/mol";
            platform::log(LogLevel::DEBUG, ss.str());
        }
    }

    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\n=== calcPairEnergy returning ===";
        ss << "\n  Final VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Final Electrostatic energy = " << elec_energy << " kJ/mol";
        ss << "\n  Total energy = " << (vdw_energy + elec_energy) << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    return {vdw_energy, elec_energy};
}

/**
 * @brief Calculate nonbonded interactions between a single residue and all other active residues
 */
inline void computeResidueNonbondedEnergy(
    model::MCState& state,
    int residue_idx,
    bool use_cutoff = false,
    bool use_pbc = false
) {
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;
    const auto& box = state.info.box;  // Box dimensions for PBC
    
    // Calculate squared cutoff distance if using cutoff
    const double cutoff2 = use_cutoff ? state.info.cutoff * state.info.cutoff : std::numeric_limits<double>::max();
    
    // Return if residue is not active
    if (!residues[residue_idx].active) {
        return;
    }
    
    // Reset energy components for current residue
    residues[residue_idx].energy_vdw = 0.0f;
    residues[residue_idx].energy_elec = 0.0f;
    
    // Iterate through all atoms in current residue
    for (int atom_i = residues[residue_idx].atomStart;
         atom_i < residues[residue_idx].atomStart + residues[residue_idx].atomCount;
         ++atom_i) {
        int type_i = atoms[atom_i].type;
        
        // Validate atom type
        if (type_i >= forcefield.numTotalTypes) {
            std::stringstream ss;
            ss << "Atom type " << type_i << " out of range. "
               << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
               << " for atom " << atom_i << " in residue " << residue_idx;
            throw std::runtime_error(ss.str());
        }
        
        // Calculate interactions with atoms in other active residues
        for (int j = 0; j < state.activeResidueCount; ++j) {
            if (!residues[j].active || j == residue_idx) continue;
            
            // Iterate through atoms in other residue
            for (int atom_j = residues[j].atomStart;
                 atom_j < residues[j].atomStart + residues[j].atomCount;
                 ++atom_j) {
                int type_j = atoms[atom_j].type;
                
                // Validate atom type
                if (type_j >= forcefield.numTotalTypes) {
                    std::stringstream ss;
                    ss << "Atom type " << type_j << " out of range. "
                       << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                       << " for atom " << atom_j << " in residue " << j;
                    throw std::runtime_error(ss.str());
                }
                
                // Calculate interatomic distance with PBC if enabled
                double dx = atoms[atom_j].x - atoms[atom_i].x;
                double dy = atoms[atom_j].y - atoms[atom_i].y;
                double dz = atoms[atom_j].z - atoms[atom_i].z;
                
                // Apply minimum image convention if PBC is enabled
                if (use_pbc) {
                    // Validate box dimensions
                    if (box[0] <= 0.0f || box[1] <= 0.0f || box[2] <= 0.0f) {
                        throw std::runtime_error("Invalid box dimensions for PBC calculation");
                    }
                    
                    // Apply minimum image convention
                    dx -= box[0] * std::round(dx / box[0]);
                    dy -= box[1] * std::round(dy / box[1]);
                    dz -= box[2] * std::round(dz / box[2]);
                    
                    if (getEnergyDebugOutput()) {
                        std::stringstream ss;
                        ss << "\nPBC distance calculation:";
                        ss << "\n  Original dx,dy,dz: " << (atoms[atom_j].x - atoms[atom_i].x)
                           << ", " << (atoms[atom_j].y - atoms[atom_i].y)
                           << ", " << (atoms[atom_j].z - atoms[atom_i].z);
                        ss << "\n  After PBC dx,dy,dz: " << dx << ", " << dy << ", " << dz;
                        ss << "\n  Box dimensions: " << box[0] << ", " << box[1] << ", " << box[2];
                        platform::log(LogLevel::DEBUG, ss.str());
                    }
                }
                
                double r2 = dx*dx + dy*dy + dz*dz;
                
                // Skip if beyond cutoff distance
                if (r2 > cutoff2) continue;
                
                // Get force field parameters
                int param_index = type_i * forcefield.numTotalTypes + type_j;
                double eps = forcefield.ljEps[param_index];
                double sigma = forcefield.ljSigma[param_index];
                double q1 = atoms[atom_i].charge;
                double q2 = atoms[atom_j].charge;
                
                // Calculate pair energy
                auto [vdw, elec] = calcPairEnergy(r2, sigma, eps, q1, q2, state.info, true);
                
                // Add energy components to current residue
                residues[residue_idx].energy_vdw += static_cast<float>(vdw);
                residues[residue_idx].energy_elec += static_cast<float>(elec);
            }
        }
    }
}

/**
 * @brief Universal function for calculating all nonbonded interactions
 */
void computeNonbondedEnergy(model::MCState& state, bool use_cutoff, bool movement_only = false, bool use_pbc = false, bool vdw_only = false) {
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

void computeSystemEnergyPBC(model::MCState& state) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation (no cutoff) ===";
        ss << "\nBox dimensions: " << state.info.box[0] << " x " 
           << state.info.box[1] << " x " << state.info.box[2] << " nm";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    // Calculate with PBC enabled but no cutoff
    computeNonbondedEnergy(state, false, false, true);
}

void computeSystemEnergyPBCCutoff(model::MCState& state) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation (with cutoff) ===";
        ss << "\nBox dimensions: " << state.info.box[0] << " x " 
           << state.info.box[1] << " x " << state.info.box[2] << " nm";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    // Calculate with both cutoff and PBC enabled
    computeNonbondedEnergy(state, true, false, true);
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

} // namespace cpu
} // namespace platform
} // namespace pygcmc 