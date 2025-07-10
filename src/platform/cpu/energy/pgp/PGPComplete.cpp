#include "PGPComplete.hpp"
#include "PGPSystem.hpp"
#include "../common/EnergyDirectCore.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include "PGPInterpolation.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Helper function declarations
static double getTotalVdwEnergy(const model::MCState& state);
static double getTotalMovementVdwEnergy(const model::MCState& state);

void computeSystemEnergyPGPComplete(model::MCState& state) {
    // First, calculate PGP electrostatic energy using existing functions
    // This includes grid interpolation, real space, and self energy
    calculateMoleculeEnergy(state);
    
    // The calculateMoleculeEnergy function already sets up ewald_energy components
    // Now we need to add complete LJ calculations
    
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
    
    // Calculate all LJ pairs with proper handling
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (!residues[i].active) continue;
        
        for (int atom_i = residues[i].atomStart;
             atom_i < residues[i].atomStart + residues[i].atomCount;
             ++atom_i) {
            
            const int type_i = atoms[atom_i].type;
            const double xi = atoms[atom_i].x;
            const double yi = atoms[atom_i].y;
            const double zi = atoms[atom_i].z;
            
            // Calculate with all other atoms (j > i to avoid double counting)
            for (int j = i; j < state.activeResidueCount; ++j) {
                if (!residues[j].active) continue;
                
                for (int atom_j = residues[j].atomStart;
                     atom_j < residues[j].atomStart + residues[j].atomCount;
                     ++atom_j) {
                    
                    // Skip self-interaction
                    if (atom_i == atom_j) continue;
                    
                    const int type_j = atoms[atom_j].type;
                    const double xj = atoms[atom_j].x;
                    const double yj = atoms[atom_j].y;
                    const double zj = atoms[atom_j].z;
                    
                    // Calculate distance with PBC
                    float dx = static_cast<float>(xi - xj);
                    float dy = static_cast<float>(yi - yj);
                    float dz = static_cast<float>(zi - zj);
                    
                    // Apply minimum image convention
                    applyPBC(dx, dy, dz, state.info.box);
                    
                    const double r2 = dx*dx + dy*dy + dz*dz;
                    
                    // Apply cutoff
                    if (r2 > cutoff2) continue;
                    
                    // Get LJ parameters
                    const int param_index = type_i * forcefield.numTotalTypes + type_j;
                    const double eps = forcefield.ljEps[param_index];
                    const double sigma = forcefield.ljSigma[param_index];
                    
                    // Calculate LJ energy
                    const double sigma2 = sigma * sigma;
                    const double sigma6 = sigma2 * sigma2 * sigma2;
                    const double sigma12 = sigma6 * sigma6;
                    const double r6 = r2 * r2 * r2;
                    const double r12 = r6 * r6;
                    
                    const double vdw = 4.0 * eps * (sigma12/r12 - sigma6/r6);
                    
                    // Add to both residues if different, or once if same residue
                    if (i == j) {
                        // Intramolecular interaction - add full energy to residue
                        residues[i].energy_vdw += static_cast<float>(vdw);
                    } else {
                        // Intermolecular - split energy between residues
                        residues[i].energy_vdw += static_cast<float>(0.5 * vdw);
                        residues[j].energy_vdw += static_cast<float>(0.5 * vdw);
                    }
                }
            }
        }
    }
    
    platform::log(LogLevel::INFO, "PGP Complete: elec=", 
                 state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self,
                 " vdw=", getTotalVdwEnergy(state));
}

void computeMovementEnergyPGPComplete(model::MCState& state) {
    // First, calculate PGP electrostatic energy for movement residues
    // using the existing movement-only functions
    
    // 1. Grid interpolation for movement residues
    double energy;
    interpolateMoleculeEnergy(state, energy);
    
    // 2. Real space for movement residues
    computeRealSpacePGP(state, true);
    
    // 3. Self energy for movement residues
    computeSelfEnergyPGP(state, true);
    
    // Now calculate LJ for movement residues
    // This needs to include interactions with ALL residues (fixed and moveable)
    computeMovementVdwEnergyDirect(state, true, true);
    
    // The energy components are already stored in state.ewald_energy
    // and residue vdw energies
    
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