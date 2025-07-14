#include "PMEReal.hpp"
#include "PMECore.hpp"
#include "platform/platform.hpp"
#include "../lj/LJMain.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

// Use COULOMB constant from energyCommon.hpp

/**
 * @brief Calculate pair energy using PME for real space
 */
std::pair<double, double> calcPairEnergyPME(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool is_excluded)
{
    double r = std::sqrt(r2);
    double lj_energy = 0.0;
    double elec_energy = 0.0;
    
    // LJ energy calculation - use unified calculateLJEnergy interface
    lj_energy = lj::calculateLJEnergyWithSwitching(r2, sigma, eps, info);
    
    // Electrostatic energy - use PME approximation
    if (!is_excluded && r < pme_params.cutoff) {
        // Normal pairs get erfc(αr)/r
        double erfc_term = pme_params.erfcApprox(r) / r;
        elec_energy = q1 * q2 * erfc_term;
    } else if (is_excluded) {
        // For excluded pairs, we need to subtract erf(αr)/r to compensate for reciprocal space
        double erfc_term = pme_params.erfcApprox(r);
        double erf_term = 1.0 - erfc_term;  // erf(x) = 1 - erfc(x)
        elec_energy = -q1 * q2 * erf_term / r;  // Note the negative sign
    }
    
    return {lj_energy, elec_energy};
}

/**
 * @brief Calculate real-space part of PME
 * 
 * @param state MC state
 * @param movement_only Whether to calculate only for moving residues
 * @param store_in_residues Whether to store energy in residues
 */
void computeRealSpacePME(model::MCState& state, bool movement_only, bool store_in_residues) {
    const auto& box = state.info.box;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float cutoff2 = pme_params.cutoff * pme_params.cutoff;

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
    const bool enable_debug = false; // Disable debug output for now

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
        
        // First handle intra-residue pairs (atoms within the same residue)
        for(int i = residues[r1].atomStart; 
            i < residues[r1].atomStart + residues[r1].atomCount - 1 && i < static_cast<int>(atoms.size()); i++) {
            if(i < 0 || i >= static_cast<int>(atoms.size()) || i >= state.activeAtomCount) continue;
            
            for(int j = i + 1; j < residues[r1].atomStart + residues[r1].atomCount && j < static_cast<int>(atoms.size()); j++) {
                if(j < 0 || j >= static_cast<int>(atoms.size()) || j >= state.activeAtomCount) continue;
                
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
                
                // Calculate real space contribution for PME - only erfc part
                double term = pme_params.erfcApprox(r);
                double pair_energy = qi * qj * term / r;
                
                // Print debug information
                if (enable_debug && debug_count < max_debug_pairs) {
                    platform::log(LogLevel::INFO, 
                        "Debug energyPME (intra): Atom pair (", i, ",", j, "): ",
                        "r = ", r, " nm, ",
                        "q1*q2 = ", qi * qj, ", ",
                        "erfc term = ", term, ", ",
                        "energy = ", pair_energy, " kJ/mol");
                    debug_count++;
                }

                // Accumulate to total energy
                real_space_total += pair_energy;
                
                // Decide how to store energy based on parameters
                if (store_in_residues) {
                    // For intra-residue pairs, all energy goes to the same residue
                    residues[r1].energy_elec += pair_energy;
                }
            }
        }
        
        // Then handle inter-residue pairs
        for(int r2 = r1 + 1; r2 < state.activeResidueCount; r2++) {
            if(!residues[r2].active) continue;
            
            // Loop over atoms in each residue with bounds checking
            for(int i = residues[r1].atomStart; 
                i < residues[r1].atomStart + residues[r1].atomCount && i < static_cast<int>(atoms.size()); i++) {
                // Ensure atom index is valid
                if(i < 0 || i >= static_cast<int>(atoms.size()) || i >= state.activeAtomCount) continue;
                
                for(int j = residues[r2].atomStart; 
                    j < residues[r2].atomStart + residues[r2].atomCount && j < static_cast<int>(atoms.size()); j++) {
                    // Ensure atom index is valid
                    if(j < 0 || j >= static_cast<int>(atoms.size()) || j >= state.activeAtomCount) continue;
                    
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
                    
                    // Calculate real space contribution for PME - only erfc part
                    double term = pme_params.erfcApprox(r);
                    double pair_energy = qi * qj * term / r;
                    
                    // Print debug information
                    if (enable_debug && debug_count < max_debug_pairs) {
                        platform::log(LogLevel::INFO, 
                            "Debug energyPME (inter): Atom pair (", i, ",", j, "): ",
                            "r = ", r, " nm, ",
                            "q1*q2 = ", qi * qj, ", ",
                            "erfc term = ", term, ", ",
                            "energy = ", pair_energy, " kJ/mol");
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
    
    // Store total real-space energy (WITHOUT COULOMB constant, to match Ewald)
    state.ewald_energy.real_space = real_space_total;
}

/**
 * @brief Compute real space PME energy for the entire system
 */
void computeRealSpaceEnergy(model::MCState& state, 
                          bool movement_only, 
                          bool store_in_residues) {
    computeRealSpacePME(state, movement_only, store_in_residues);
}

/**
 * @brief Calculate interaction energy between two specific atoms using PME real space
 */
double calculateAtomPairEnergyRealSpace(const model::MCAtom& atom1, const model::MCAtom& atom2,
                                      const float box[3], double cutoff) {
    // Calculate distance
    float dx = atom1.x - atom2.x;
    float dy = atom1.y - atom2.y;
    float dz = atom1.z - atom2.z;
    
    // Apply periodic boundary conditions
    dx -= box[0] * round(dx / box[0]);
    dy -= box[1] * round(dy / box[1]);
    dz -= box[2] * round(dz / box[2]);
    
    float r2 = dx*dx + dy*dy + dz*dz;
    float r = sqrt(r2);
    
    // Return 0 if beyond cutoff
    if (r >= cutoff) return 0.0;
    
    // Skip neutral atoms
    if (std::abs(atom1.charge) < 1e-6 || std::abs(atom2.charge) < 1e-6) return 0.0;
    
    // Calculate real space contribution - only erfc part
    double term = pme_params.erfcApprox(r);
    double energy = atom1.charge * atom2.charge * term / r;
    
    return energy;
}

/**
 * @brief Validate real space parameters
 */
bool validateRealSpaceParameters(double cutoff, double alpha) {
    if (cutoff <= 0.0) {
        platform::log(LogLevel::ERROR, "Invalid cutoff: ", cutoff);
        return false;
    }
    
    if (alpha <= 0.0) {
        platform::log(LogLevel::ERROR, "Invalid alpha: ", alpha);
        return false;
    }
    
    // Check if alpha is reasonable for the given cutoff
    double alphaR = alpha * cutoff;
    if (alphaR < 2.0) {
        platform::log(LogLevel::WARNING, "Alpha may be too small for cutoff. alphaR = ", alphaR);
    }
    
    if (alphaR > 10.0) {
        platform::log(LogLevel::WARNING, "Alpha may be too large for cutoff. alphaR = ", alphaR);
    }
    
    return true;
}

// <agent-hook:real_space_implementation>

} // namespace cpu
} // namespace platform
} // namespace pygcmc