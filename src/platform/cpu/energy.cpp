// src/platform/cpu/energy.cpp
#include "energy.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>  // For output formatting

namespace pygcmc {
namespace platform {
namespace cpu {

// Debug flag to control output - default to false for production use
static bool debug_output = false;

/**
 * @brief Coulomb constant in GROMACS MD units [kJ·nm/mol/e²]
 * 
 * k_c = 1/(4*π*ε₀) = 138.935456 kJ·nm/mol/e²
 * 
 * Unit analysis:
 * - ε₀ (vacuum permittivity) = 8.8541878128e-12 C²/(J·m)
 * - 1 kJ = 1000 J
 * - 1 nm = 1e-9 m
 * - 1 e = 1.60217663e-19 C
 * - N_A (Avogadro constant) = 6.02214076e23 mol⁻¹
 */
const float COULOMB = 138.935456f;

/**
 * @brief Safety parameters for energy calculation
 * 
 * !!! CRITICAL: Distance handling for Monte Carlo simulation !!!
 * 
 * MIN_SAFE_DISTANCE: Minimum allowed distance (1% of sigma)
 * - !!! Prevents numerical instability and infinity at r = 0
 * - !!! Essential for Monte Carlo sampling near contact
 * - !!! Implements soft core potential for r < MIN_SAFE_DISTANCE
 * 
 * MAX_SAFE_ENERGY: Maximum allowed energy per interaction
 * - !!! Prevents numerical overflow in Metropolis criterion
 * - !!! Keeps energies finite for stable MC sampling
 * - !!! Especially important for Coulomb interactions at small r
 */
const float MIN_SAFE_DISTANCE = 0.01f;  // nm (1% of typical sigma)
const float MAX_SAFE_ENERGY = 1e6f;     // kJ/mol

/**
 * @brief Calculate LJ and Coulomb energy with safety checks
 * 
 * !!! IMPORTANT: Zero distance handling strategy !!!
 * 1. For r < MIN_SAFE_DISTANCE:
 *    - Replace actual distance with MIN_SAFE_DISTANCE
 *    - Provides continuous potential without singularity
 *    - Allows MC moves through high-energy regions
 * 
 * 2. Energy capping:
 *    - Limits maximum repulsion to MAX_SAFE_ENERGY
 *    - Prevents exp(−βE) underflow in Metropolis
 *    - Maintains numerical stability of MC sampling
 * 
 * This approach:
 * - !!! Avoids infinite energies at r = 0
 * - !!! Keeps energy continuous and differentiable
 * - !!! Allows MC sampling of close contacts
 * - !!! Prevents numerical instabilities in simulation
 */
inline std::pair<float, float> calcPairEnergy(float r2, float sigma, float eps, float q1, float q2) {
    if (debug_output) {
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

    // Apply minimum safe distance for numerical stability
    if (r2 < MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE) {
        if (debug_output) {
            platform::log(LogLevel::DEBUG, "Distance below MIN_SAFE_DISTANCE, using r2 = ", 
                         MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE);
        }
        r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    }
    
    float r = std::sqrt(r2);
    
    if (debug_output) {
        platform::log(LogLevel::DEBUG, "Distance r = ", r, " nm");
    }
    
    // Calculate LJ energy: V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    // 使用更稳定的计算方法：先计算(σ/r)²，然后通过乘法得到6次和12次方
    float sigma_r2 = (sigma * sigma) / r2;  // (σ/r)²
    float sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;  // (σ/r)⁶
    float sigma_r12 = sigma_r6 * sigma_r6;  // (σ/r)¹²
    float vdw_energy = 4.0f * eps * (sigma_r12 - sigma_r6);  // kJ/mol
    
    // Calculate Coulomb energy: V_C = k_c * q1*q2/r
    float elec_energy = COULOMB * q1 * q2 / r;  // kJ/mol

    if (debug_output) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\nEnergy calculation details:";
        ss << "\n  (sigma/r)² = " << sigma_r2;
        ss << "\n  (sigma/r)⁶ = " << sigma_r6;
        ss << "\n  (sigma/r)¹² = " << sigma_r12;
        ss << "\n  4*epsilon = " << (4.0f * eps);
        ss << "\n  VDW term = " << (sigma_r12 - sigma_r6);
        ss << "\n  COULOMB constant = " << COULOMB;
        ss << "\n  q1*q2 = " << (q1 * q2);
        ss << "\nInitial energies:";
        ss << "\n  VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Electrostatic energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // !!! CRITICAL: Apply energy capping for numerical stability
    // First cap individual terms
    if (debug_output && (std::abs(vdw_energy) > MAX_SAFE_ENERGY || std::abs(elec_energy) > MAX_SAFE_ENERGY)) {
        std::stringstream ss;
        ss << "\nEnergy capping applied:";
        ss << "\n  Original VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Original Elec energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    vdw_energy = std::min(vdw_energy, MAX_SAFE_ENERGY);
    vdw_energy = std::max(vdw_energy, -MAX_SAFE_ENERGY);
    elec_energy = std::min(elec_energy, MAX_SAFE_ENERGY);
    elec_energy = std::max(elec_energy, -MAX_SAFE_ENERGY);
    
    if (debug_output && (std::abs(vdw_energy) > MAX_SAFE_ENERGY || std::abs(elec_energy) > MAX_SAFE_ENERGY)) {
        std::stringstream ss;
        ss << "\nAfter individual capping:";
        ss << "\n  Capped VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Capped Elec energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // !!! CRITICAL: Also cap total energy
    float total_energy = vdw_energy + elec_energy;
    float original_total = total_energy;
    if (total_energy > MAX_SAFE_ENERGY) {
        float scale = MAX_SAFE_ENERGY / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (debug_output) {
            std::stringstream ss;
            ss << "\nTotal energy exceeded MAX_SAFE_ENERGY:";
            ss << "\n  Original total = " << original_total << " kJ/mol";
            ss << "\n  Scale factor = " << scale;
            ss << "\n  Final VDW = " << vdw_energy << " kJ/mol";
            ss << "\n  Final Elec = " << elec_energy << " kJ/mol";
            ss << "\n  Final total = " << (vdw_energy + elec_energy) << " kJ/mol";
            platform::log(LogLevel::DEBUG, ss.str());
        }
    } else if (total_energy < -MAX_SAFE_ENERGY) {
        float scale = -MAX_SAFE_ENERGY / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (debug_output) {
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

    if (debug_output) {
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
 * 
 * @param state System state
 * @param residue_idx Index of the residue to calculate energy for
 * @param use_cutoff Whether to use distance cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 * 
 * This function handles three scenarios:
 * 1. No cutoff, no PBC: Calculate interactions between all atom pairs
 * 2. With cutoff, no PBC: Only calculate interactions within cutoff distance
 * 3. With PBC: Apply minimum image convention for distance calculation
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
    const float cutoff2 = use_cutoff ? state.info.cutoff * state.info.cutoff : std::numeric_limits<float>::max();
    
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
                float dx = atoms[atom_j].x - atoms[atom_i].x;
                float dy = atoms[atom_j].y - atoms[atom_i].y;
                float dz = atoms[atom_j].z - atoms[atom_i].z;
                
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
                    
                    if (debug_output) {
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
                
                float r2 = dx*dx + dy*dy + dz*dz;
                
                // Skip if beyond cutoff distance
                if (r2 > cutoff2) continue;
                
                // Get force field parameters
                int param_index = type_i * forcefield.numTotalTypes + type_j;
                float eps = forcefield.ljEps[param_index];
                float sigma = forcefield.ljSigma[param_index];
                float q1 = atoms[atom_i].charge;
                float q2 = atoms[atom_j].charge;
                
                // Calculate pair energy
                auto [vdw, elec] = calcPairEnergy(r2, sigma, eps, q1, q2);
                
                // Add energy components to current residue
                residues[residue_idx].energy_vdw += vdw;
                residues[residue_idx].energy_elec += elec;
            }
        }
    }
}

/**
 * @brief Universal function for calculating all nonbonded interactions
 * 
 * @param state System state
 * @param use_cutoff Whether to use distance cutoff
 * @param movement_only Whether to calculate only for movement residues
 * @param use_pbc Whether to use periodic boundary conditions
 * 
 * This function provides a unified interface for all nonbonded energy calculations:
 * - Can handle both cutoff and non-cutoff calculations
 * - Can calculate for all residues or movement residues only
 * - Validates force field parameters based on calculation type
 * - Provides detailed debug output for energy components
 */
void computeNonbondedEnergy(model::MCState& state, bool use_cutoff, bool movement_only = false, bool use_pbc = false) {
    if (debug_output) {
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
        residue.energy_elec = 0.0f;
    }

    if (movement_only) {
        // Validate movement residues
        if (state.movementResidues.empty()) {
            throw std::runtime_error("No movement residues defined");
        }

        // Calculate energies only for movement residues
        for (const auto& movementInfo : state.movementResidues) {
            if (debug_output) {
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

                if (debug_output) {
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
    if (debug_output) {
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

/**
 * @brief Interface functions for nonbonded energy calculations
 */

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

/**
 * @brief Interface functions for nonbonded energy calculations with PBC support
 * 
 * These functions provide PBC-specific versions of the energy calculations:
 * - Uses minimum image convention for distance calculations
 * - Requires valid box dimensions in state.info.box
 * - Recommended to use with cutoff for better performance
 */

void computeSystemEnergyPBC(model::MCState& state) {
    if (debug_output) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation ===";
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

// Function to enable/disable debug output
void setEnergyDebugOutput(bool enable) {
    debug_output = enable;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc


