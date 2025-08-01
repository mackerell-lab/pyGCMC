#include "PGPComplete.hpp"
#include "PGPGlobal.hpp"
#include "PGPSystem.hpp"
#include "PGPCore.hpp"
#include "PGPReal.hpp"
#include "PGPSelf.hpp"
#include "PGPInterpolation.hpp"
#include "../common/EnergyDirectCore.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Helper function declarations
static void calculateCompleteLJEnergy(model::MCState& state);
static void calculateCompleteRealSpaceElectrostatics(model::MCState& state);
static double getTotalVdwEnergy(const model::MCState& state);
static double getTotalMovementVdwEnergy(const model::MCState& state);

void computeSystemEnergyPGPComplete(model::MCState& state) {
    // PGP Complete uses pure PGP method for electrostatics
    // but includes ALL interactions (including intramolecular)
    
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    platform::log(LogLevel::DEBUG, "Computing complete system energy using PGP method");
    
    // 1. Calculate grid potential interpolation part (PGP interpolation)
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate complete real space electrostatics (including intramolecular)
    calculateCompleteRealSpaceElectrostatics(state);
    
    // 3. Calculate self energy correction
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, false);
    
    // 4. Calculate complete LJ interactions (including intramolecular)
    calculateCompleteLJEnergy(state);
    
    // Apply COULOMB constant to real space energy
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate total energies
    double vdw_total = getTotalVdwEnergy(state);
    
    // Store grid energy in reciprocal field for consistency
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    platform::log(LogLevel::INFO, "PGP Complete: grid=", grid_energy,
                 " real=", state.ewald_energy.real_space,
                 " self=", state.ewald_energy.self,
                 " vdw=", vdw_total);
}

void computeMovementEnergyPGPComplete(model::MCState& state) {
    // For movement residues, calculate complete interactions
    
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    platform::log(LogLevel::DEBUG, "Computing movement energy using PGP Complete method");
    
    // 1. Calculate grid potential interpolation for movement residues
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. Calculate real space electrostatics for movement residues
    // This needs to include ALL interactions of movement atoms
    calculateCompleteRealSpaceElectrostatics(state);
    
    // 3. Calculate self energy for movement residues
    state.ewald_energy.self = computeSelfEnergyPGPImpl(state, true);
    
    // 4. Calculate complete LJ for movement residues
    // Reset VDW energies for movement residues first
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                state.residues[i].energy_vdw = 0.0f;
            }
        }
    }
    
    // Calculate LJ interactions for movement residues with all atoms
    calculateCompleteLJEnergy(state);
    
    // Apply COULOMB constant
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate totals
    double vdw_total = getTotalMovementVdwEnergy(state);
    
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    platform::log(LogLevel::INFO, "PGP Complete Movement: grid=", grid_energy,
                 " real=", state.ewald_energy.real_space,  
                 " self=", state.ewald_energy.self,
                 " vdw=", vdw_total);
}

// Calculate complete real space electrostatics including intramolecular
static void calculateCompleteRealSpaceElectrostatics(model::MCState& state) {
    const auto& atoms = state.atoms;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    const double alpha = getPGPParams().alpha;
    
    // Reset real space energy
    state.ewald_energy.real_space = 0.0;
    
    // Calculate all pairwise electrostatic interactions
    for (int i = 0; i < state.activeAtomCount - 1; ++i) {
        const double qi = atoms[i].charge;
        if (qi == 0.0) continue;
        
        const double xi = atoms[i].x;
        const double yi = atoms[i].y;
        const double zi = atoms[i].z;
        
        for (int j = i + 1; j < state.activeAtomCount; ++j) {
            const double qj = atoms[j].charge;
            if (qj == 0.0) continue;
            
            // Calculate distance with PBC
            double dx = xi - atoms[j].x;
            double dy = yi - atoms[j].y;
            double dz = zi - atoms[j].z;
            
            // Apply minimum image convention
            float dx_f = static_cast<float>(dx);
            float dy_f = static_cast<float>(dy);
            float dz_f = static_cast<float>(dz);
            applyPBC(dx_f, dy_f, dz_f, state.info.box);
            dx = dx_f; dy = dy_f; dz = dz_f;
            
            const double r2 = dx*dx + dy*dy + dz*dz;
            
            // Apply cutoff
            if (r2 > cutoff2) continue;
            
            // Skip extremely close atoms
            if (r2 < 1e-12) continue;
            
            const double r = std::sqrt(r2);
            
            // Calculate erfc(alpha*r)/r
            const double alphar = alpha * r;
            const double erfc_val = std::erfc(alphar);
            
            // Real space electrostatic energy
            const double energy = qi * qj * erfc_val / r;
            
            // Skip if not finite
            if (!std::isfinite(energy)) continue;
            
            state.ewald_energy.real_space += energy;
        }
    }
}

// Calculate complete LJ energy including intramolecular interactions
static void calculateCompleteLJEnergy(model::MCState& state) {
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;
    const double cutoff2 = state.info.cutoff * state.info.cutoff;
    
    // Reset all VDW energies
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
    }
    
    // Create atom to residue mapping for efficiency
    std::vector<int> atomToResidue(state.activeAtomCount, -1);
    for (int r = 0; r < state.activeResidueCount; ++r) {
        if (!residues[r].active) continue;
        for (int a = residues[r].atomStart; 
             a < residues[r].atomStart + residues[r].atomCount; ++a) {
            atomToResidue[a] = r;
        }
    }
    
    // Calculate all LJ pairs without double counting
    for (int atom_i = 0; atom_i < state.activeAtomCount - 1; ++atom_i) {
        int res_i = atomToResidue[atom_i];
        if (res_i == -1) continue;
        
        const int type_i = atoms[atom_i].type;
        const double xi = atoms[atom_i].x;
        const double yi = atoms[atom_i].y;
        const double zi = atoms[atom_i].z;
        
        for (int atom_j = atom_i + 1; atom_j < state.activeAtomCount; ++atom_j) {
            int res_j = atomToResidue[atom_j];
            if (res_j == -1) continue;
            
            const int type_j = atoms[atom_j].type;
            
            // Calculate distance with PBC
            double dx = xi - atoms[atom_j].x;
            double dy = yi - atoms[atom_j].y;
            double dz = zi - atoms[atom_j].z;
            
            // Apply minimum image convention
            float dx_f = static_cast<float>(dx);
            float dy_f = static_cast<float>(dy);
            float dz_f = static_cast<float>(dz);
            applyPBC(dx_f, dy_f, dz_f, state.info.box);
            dx = dx_f; dy = dy_f; dz = dz_f;
            
            const double r2 = dx*dx + dy*dy + dz*dz;
            
            // Apply cutoff
            if (r2 > cutoff2) continue;
            
            // Skip extremely close atoms
            if (r2 < 1e-12) continue;
            
            // Get LJ parameters
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
            
            // Skip if not finite
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