#include "PGPComplete.hpp"
#include "PGPSystem.hpp"
#include "../common/EnergyDirectCore.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include "PGPInterpolation.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "../pme/PMEComposite.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Helper function declarations
static double getTotalVdwEnergy(const model::MCState& state);
static double getTotalMovementVdwEnergy(const model::MCState& state);

void computeSystemEnergyPGPComplete(model::MCState& state) {
    // For PGP Complete, we use PME for electrostatics (which already includes all interactions)
    // and only recalculate LJ to include intramolecular interactions
    
    // First calculate standard PME electrostatics
    platform::cpu::PMEComposite::computeSystemEnergy(state);
    
    // Now we need to recalculate LJ to include intramolecular interactions
    
    // Calculate ALL LJ interactions, including intramolecular
    // We'll use a modified approach that includes all pairs
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    
    // Reset VDW energies
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
    }
    
    // Calculate all LJ pairs without double counting
    // Use atom-based loops to ensure each pair is counted exactly once
    for (int atom_i = 0; atom_i < state.activeAtomCount - 1; ++atom_i) {
        // Find which residue atom_i belongs to
        int res_i = -1;
        for (int r = 0; r < state.activeResidueCount; ++r) {
            if (!residues[r].active) continue;
            if (atom_i >= residues[r].atomStart && 
                atom_i < residues[r].atomStart + residues[r].atomCount) {
                res_i = r;
                break;
            }
        }
        if (res_i == -1) continue;  // Skip if atom not in active residue
        
        const int type_i = atoms[atom_i].type;
        const double xi = atoms[atom_i].x;
        const double yi = atoms[atom_i].y;
        const double zi = atoms[atom_i].z;
        
        // Only check atoms after atom_i to avoid double counting
        for (int atom_j = atom_i + 1; atom_j < state.activeAtomCount; ++atom_j) {
            // Find which residue atom_j belongs to
            int res_j = -1;
            for (int r = 0; r < state.activeResidueCount; ++r) {
                if (!residues[r].active) continue;
                if (atom_j >= residues[r].atomStart && 
                    atom_j < residues[r].atomStart + residues[r].atomCount) {
                    res_j = r;
                    break;
                }
            }
            if (res_j == -1) continue;  // Skip if atom not in active residue
            
            const int type_j = atoms[atom_j].type;
            const double xj = atoms[atom_j].x;
            const double yj = atoms[atom_j].y;
            const double zj = atoms[atom_j].z;
            
            // Calculate distance with PBC
            double dx = xi - xj;
            double dy = yi - yj;
            double dz = zi - zj;
            
            // Apply minimum image convention
            float dx_f = static_cast<float>(dx);
            float dy_f = static_cast<float>(dy);
            float dz_f = static_cast<float>(dz);
            applyPBC(dx_f, dy_f, dz_f, state.info.box);
            dx = dx_f; dy = dy_f; dz = dz_f;
            
            const double r2 = dx*dx + dy*dy + dz*dz;
            
            // Apply cutoff
            if (r2 > cutoff2) continue;
            
            // Skip extremely close atoms to avoid divide-by-zero
            if (r2 < 1e-12) continue;
            
            // Bounds check for LJ parameters
            const int param_index = type_i * forcefield.numTotalTypes + type_j;
            if (param_index < 0 || param_index >= static_cast<int>(forcefield.ljEps.size())) {
                platform::log(LogLevel::WARNING, "Invalid LJ parameter index: ", param_index);
                continue;
            }
            
            const double eps = forcefield.ljEps[param_index];
            const double sigma = forcefield.ljSigma[param_index];
            
            // Skip if no LJ interaction
            if (eps == 0.0 || sigma == 0.0) continue;
            
            // Calculate LJ energy
            const double sigma2 = sigma * sigma;
            const double sigma6 = sigma2 * sigma2 * sigma2;
            const double sigma12 = sigma6 * sigma6;
            const double r6 = r2 * r2 * r2;
            const double r12 = r6 * r6;
            
            const double vdw = 4.0 * eps * (sigma12/r12 - sigma6/r6);
            
            // Skip if result is not finite (NaN or Inf)
            if (!std::isfinite(vdw)) continue;
            
            // Add energy to residues
            if (res_i == res_j) {
                // Intramolecular interaction - add full energy to the residue
                residues[res_i].energy_vdw += static_cast<float>(vdw);
            } else {
                // Intermolecular - split energy between residues
                residues[res_i].energy_vdw += static_cast<float>(0.5 * vdw);
                residues[res_j].energy_vdw += static_cast<float>(0.5 * vdw);
            }
        }
    }
    
    platform::log(LogLevel::INFO, "PGP Complete: elec=", 
                 state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self,
                 " vdw=", getTotalVdwEnergy(state));
}

void computeMovementEnergyPGPComplete(model::MCState& state) {
    // For PGP Complete movement energy, we use PME for electrostatics
    // (since PGP interpolation doesn't capture all interactions properly)
    // and calculate complete LJ including intramolecular
    
    // First calculate PME movement energy for electrostatics
    platform::cpu::PMEComposite::computeMovementEnergy(state);
    
    // Now calculate LJ for movement residues
    // Reset VDW energies for movement residues
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                state.residues[i].energy_vdw = 0.0f;
            }
        }
    }
    
    // Calculate LJ interactions for movement residues with all atoms
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    
    // For each movement residue
    for (const auto& movementInfo : state.movementResidues) {
        for (int res_idx = movementInfo.startIndex;
             res_idx < movementInfo.startIndex + movementInfo.activeCount; res_idx++) {
            if (!residues[res_idx].active) continue;
            
            // For each atom in this movement residue
            for (int atom_i = residues[res_idx].atomStart;
                 atom_i < residues[res_idx].atomStart + residues[res_idx].atomCount;
                 ++atom_i) {
                
                const int type_i = atoms[atom_i].type;
                const double xi = atoms[atom_i].x;
                const double yi = atoms[atom_i].y;
                const double zi = atoms[atom_i].z;
                
                // Check against ALL other atoms (including fixed residues)
                for (int atom_j = 0; atom_j < state.activeAtomCount; ++atom_j) {
                    if (atom_i == atom_j) continue;
                    
                    // Find which residue atom_j belongs to
                    int res_j = -1;
                    for (int r = 0; r < state.activeResidueCount; ++r) {
                        if (!residues[r].active) continue;
                        if (atom_j >= residues[r].atomStart && 
                            atom_j < residues[r].atomStart + residues[r].atomCount) {
                            res_j = r;
                            break;
                        }
                    }
                    if (res_j == -1) continue;
                    
                    const int type_j = atoms[atom_j].type;
                    const double xj = atoms[atom_j].x;
                    const double yj = atoms[atom_j].y;
                    const double zj = atoms[atom_j].z;
                    
                    // Calculate distance with PBC
                    double dx = xi - xj;
                    double dy = yi - yj;
                    double dz = zi - zj;
                    
                    // Apply minimum image convention
                    float dx_f = static_cast<float>(dx);
                    float dy_f = static_cast<float>(dy);
                    float dz_f = static_cast<float>(dz);
                    applyPBC(dx_f, dy_f, dz_f, state.info.box);
                    dx = dx_f; dy = dy_f; dz = dz_f;
                    
                    const double r2 = dx*dx + dy*dy + dz*dz;
                    
                    // Apply cutoff
                    if (r2 > cutoff2) continue;
                    
                    // Skip extremely close atoms to avoid divide-by-zero
                    if (r2 < 1e-12) continue;
                    
                    // Get LJ parameters
                    const int param_index = type_i * forcefield.numTotalTypes + type_j;
                    if (param_index < 0 || param_index >= static_cast<int>(forcefield.ljEps.size())) {
                        continue;
                    }
                    
                    const double eps = forcefield.ljEps[param_index];
                    const double sigma = forcefield.ljSigma[param_index];
                    
                    if (eps == 0.0 || sigma == 0.0) continue;
                    
                    // Calculate LJ energy
                    const double sigma2 = sigma * sigma;
                    const double sigma6 = sigma2 * sigma2 * sigma2;
                    const double sigma12 = sigma6 * sigma6;
                    const double r6 = r2 * r2 * r2;
                    const double r12 = r6 * r6;
                    
                    const double vdw = 4.0 * eps * (sigma12/r12 - sigma6/r6);
                    
                    // Skip if result is not finite (NaN or Inf)
                    if (!std::isfinite(vdw)) continue;
                    
                    // For movement energy, add full energy to movement residue
                    // (the other residue's contribution will be calculated when it moves)
                    residues[res_idx].energy_vdw += static_cast<float>(vdw);
                }
            }
        }
    }
    
    platform::log(LogLevel::INFO, "PGP Complete Movement: elec=", 
                 state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self,
                 " vdw=", getTotalMovementVdwEnergy(state));
}

// Helper function to get total VDW energy
static double getTotalVdwEnergy(const model::MCState& state) {
    double total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            total += residue.energy_vdw;
        }
    }
    return total;
}

// Helper function to get movement VDW energy
static double getTotalMovementVdwEnergy(const model::MCState& state) {
    double total = 0.0;
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                total += state.residues[i].energy_vdw;
            }
        }
    }
    return total;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc