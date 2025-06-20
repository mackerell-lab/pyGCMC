#include "PGPEvaluator.hpp"
#include <cmath>
#include <iostream>

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
    const float cutoff2 = pgp_params.cutoff * pgp_params.cutoff;

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
            
            // 如果两个残基都是固定的，那么它们之间的相互作用已经包含在预计算的网格势能中
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
                    double term = pgp_params.erfcApprox(r);
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
 * @brief Calculate self energy correction for PGP method
 */
double computeSelfEnergyPGPImpl(model::MCState& state, bool movement_only) {
    double self_energy = 0.0;
    double sum_q2 = 0.0;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing self energy for PGP method");
        platform::log(LogLevel::DEBUG, "Movement only: ", movement_only);
    }
    
    // Sum up squares of charges
    if (movement_only) {
        // Only include moving residues
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex; 
                 i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if (!state.residues[i].active) continue;
                
                // Sum q^2 for all atoms in this movement residue
                for (int j = 0; j < state.residues[i].atomCount; j++) {
                    int atomIdx = state.residues[i].atomStart + j;
                    double q = state.atoms[atomIdx].charge;
                    sum_q2 += q * q;
                }
            }
        }
    } 
    else {
        // Sum q^2 for all atoms in the system
        for(int i = 0; i < state.activeAtomCount; i++) {
            double charge = state.atoms[i].charge;
            double q2 = charge * charge;
            sum_q2 += q2;
        }
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Sum of q² = ", sum_q2);
    }
    
    // Self-energy formula: -ONE_4PI_EPS0 * alpha / sqrt(PI) * sum_q2
    double prefactor = -COULOMB * pgp_params.alpha / sqrt(M_PI);
    self_energy = prefactor * sum_q2;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Self energy prefactor = ", prefactor, 
                     ", resulting self energy = ", self_energy);
    }
    
    return self_energy;
}

/**
 * @brief Use PGP method to calculate system energy
 */
void computeSystemEnergyPGPImpl(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing total system energy using PGP method");
    }
    
    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space part
    computeRealSpacePGPImpl(state, false, true);
    
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
        platform::log(LogLevel::DEBUG, "PGP system energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

/**
 * @brief Use PGP method to calculate energy of moving residues
 */
void computeMovementEnergyPGPImpl(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing movement residue energy using PGP method");
    }
    
    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space part (only for moving residues)
    computeRealSpacePGPImpl(state, true, true);
    
    // 3. Calculate self energy correction (only for moving residues)
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, true);
    
    // 4. Calculate LJ interactions using direct cutoff method
    computeSystemVdwEnergyCutoff(state);
    
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
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP movement energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

// Public interface wrappers
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues) {
    computeRealSpacePGPImpl(state, movement_only, store_in_residues);
}

double computeSelfEnergyPGP(model::MCState& state, bool movement_only) {
    return computeSelfEnergyPGPImpl(state, movement_only);
}

void computeSystemEnergyPGP(model::MCState& state) {
    computeSystemEnergyPGPImpl(state);
}

void computeMovementEnergyPGP(model::MCState& state) {
    computeMovementEnergyPGPImpl(state);
}

// <agent-hook:pgp_evaluator_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 