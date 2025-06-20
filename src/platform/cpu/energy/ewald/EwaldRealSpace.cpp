#include "EwaldRealSpace.hpp"
#include "../lj/LJSwitching.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate pair energy for Ewald real-space part
 */
std::pair<double, double> calcPairEnergyEwaldRealSpace(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool is_excluded)
{    
    // Calculate VdW energy using simplified interface
    double vdw_energy = lj::calcLJEnergyWithSwitching(r2, sigma, eps, info);
    
    // Apply minimum safe distance
    if (r2 < MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE) {
        r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    }
    
    double r = std::sqrt(r2);
    
    // Calculate electrostatic energy
    double elec_energy;
    if (is_excluded) {
        // For excluded pairs, subtract erf(αr)/r to compensate for reciprocal space
        double erfc_term = ewald_params.erfcApprox(r);
        double erf_term = 1.0 - erfc_term;  // erf(x) = 1 - erfc(x)
        elec_energy = -COULOMB * q1 * q2 * erf_term / r;  // Note the negative sign
    } else {
        // Normal pairs get erfc(αr)/r
        double erfc_term = ewald_params.erfcApprox(r) / r;
        elec_energy = q1 * q2 * erfc_term;
    }
    
    // Apply energy limits
    const double max_safe_energy = static_cast<double>(MAX_SAFE_ENERGY);
    elec_energy = std::min(std::max(elec_energy, -max_safe_energy), max_safe_energy);
    
    return {vdw_energy, elec_energy};
}

/**
 * @brief Calculate real-space part of Ewald summation
 */
void computeRealSpaceEwald(model::MCState& state, bool movement_only, bool store_in_residues) {
    const auto& box = state.info.box;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float cutoff2 = ewald_params.cutoff * ewald_params.cutoff;

    // Reset electrostatic energy
    for(auto& residue : residues) {
        if(residue.active) {
            residue.energy_elec = 0.0f;
        }
    }
    
    // Real-space total energy
    double real_space_total = 0.0;
    
    // Debug information
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

            // Loop over atom pairs between residues
            for(int i = residues[r1].atomStart; 
                i < residues[r1].atomStart + residues[r1].atomCount; i++) {
                
                for(int j = residues[r2].atomStart;
                    j < residues[r2].atomStart + residues[r2].atomCount; j++) {
                    
                    // Calculate minimum image distance
                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;

                    // Apply periodic boundary conditions
                    if(dx > box[0]/2) dx -= box[0];
                    else if(dx < -box[0]/2) dx += box[0];
                    if(dy > box[1]/2) dy -= box[1];
                    else if(dy < -box[1]/2) dy += box[1];
                    if(dz > box[2]/2) dz -= box[2];
                    else if(dz < -box[2]/2) dz += box[2];

                    float r2 = dx*dx + dy*dy + dz*dz;

                    // Only calculate for pairs within cutoff range
                    if(r2 < cutoff2) {
                        float r = std::sqrt(r2);
                        float qi = atoms[i].charge;
                        float qj = atoms[j].charge;

                        // Calculate erfc(αr)/r directly
                        double alphaR = ewald_params.alpha * r;
                        double term = std::erfc(alphaR) / r;
                        
                        // Calculate energy contribution
                        double pair_energy = qi * qj * term;

                        // Debug output for first few pairs
                        if (debug_count < max_debug_pairs) {
                            platform::log(LogLevel::INFO, 
                                "Debug EwaldRealSpace: Atom pair (", i, ",", j, "): ",
                                "r = ", r, " nm, ",
                                "q1*q2 = ", qi * qj, ", ",
                                "erfc term = ", term, ", ",
                                "energy = ", pair_energy,
                                ", with COULOMB = ", COULOMB * pair_energy, " kJ/mol");
                            debug_count++;
                        }

                        // Accumulate to total energy
                        real_space_total += pair_energy;
                        
                        // Store energy in residues if requested
                        if (store_in_residues) {
                            // Each residue gets half of the pair interaction energy
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
 * @brief Validate real space calculation parameters
 */
bool validateRealSpaceParameters(const model::MCState& state) {
    // Check if Ewald parameters are initialized
    if (!ewald_params.initialized) {
        platform::log(LogLevel::ERROR, "Ewald parameters not initialized");
        return false;
    }
    
    // Check periodic boundary conditions
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        platform::log(LogLevel::ERROR, "Ewald method requires periodic boundary conditions");
        return false;
    }
    
    // Check cutoff vs box size
    float minBoxSize = std::min(state.info.box[0], std::min(state.info.box[1], state.info.box[2]));
    if (ewald_params.cutoff >= 0.5f * minBoxSize) {
        platform::log(LogLevel::WARNING, 
            "Warning: Cutoff distance (", ewald_params.cutoff, 
            " nm) is larger than half the smallest box dimension (", 
            minBoxSize/2, " nm). This may affect minimum image convention.");
    }
    
    return true;
}

/**
 * @brief Get real space energy breakdown by residue
 */
void getRealSpaceEnergyBreakdown(const model::MCState& state, 
                               std::vector<double>& energies) {
    energies.clear();
    energies.reserve(state.activeResidueCount);
    
    for(int i = 0; i < state.activeResidueCount; i++) {
        if(state.residues[i].active) {
            energies.push_back(static_cast<double>(state.residues[i].energy_elec));
        } else {
            energies.push_back(0.0);
        }
    }
}

// <agent-hook:ewald_real_space_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 